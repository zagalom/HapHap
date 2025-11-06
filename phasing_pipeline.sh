#!/bin/bash

# Arquivo: phasing_and_haplotype_pipeline.sh
# Propósito: Automatizar mapeamento, faseamento (GATK + WhatsHap) e
# reconstrução de haplótipos (Python) usando ficheiros FASTQ/BAM.

set -euo pipefail

# --- CORREÇÃO CRÍTICA DE LOCALE: Garante que 'printf' e 'bc' usam o ponto (.) como decimal ---
export LC_NUMERIC="C"

# ---- CONFIGURATION ----

# Permite que o utilizador forneça entrada ou solicita se não for fornecida (com valores padrão)
GENES_FASTA=${1:-}
RESULTS_DIR=${2:-}
FASTQ_DIR=${3:-}
THREADS=${4:-4}      # Padrão: 4 threads
GENES_LIST="${5:-}"  # Argumento opcional, e.g. ERG11,GSC1,UPC2

# FATOR DE MULTIPLICAÇÃO PARA O FILTRO DP MÁXIMO (DP_MÉDIA * FATOR)
# Usamos 2X para atenuar o gene collapse e remover apenas os outliers extremos.
DP_MULTIPLIER="2"

# ATENÇÃO: Defina esta variável como "true" para IGNORAR a filtragem de Profundidade Alélica (AD).
BYPASS_FILTER="false"

# Variáveis para Read Group (RG) - Inseridas no cabeçalho pelo BWA
RG_LB="LIB1"
RG_PL="ILLUMINA"
RG_PU="L001"

HAPLOTYPE_SCRIPT="haplotype_reconstruction.py"
TEMP_VCF_SUFFIX=".unphased.vcf.gz"
FILTERED_VCF_SUFFIX=".filtered.vcf.gz" # VCF temporário filtrado (PASS)
FINAL_VCF_SUFFIX=".phased.vcf.gz"
IMBALANCE_REPORT_SUFFIX=".imbalanced_report.txt"
ALL_VARIANTS_REPORT_SUFFIX=".all_variants_filter_report.txt"

function print_usage {
    echo "Uso:"
    echo "  $0 <genes_fasta> <results_dir> <fastq_dir> [threads] [genes_list]"
    echo "    <genes_fasta>  : Ficheiro FASTA de genes de referência"
    echo "    <results_dir>  : Diretório para output/subdiretórios de amostras"
    echo "    <fastq_dir>    : Diretório contendo FASTQs emparelhados (SAMPLE_1.fastq.gz/SAMPLE_2.fastq.gz)"
    echo "    [threads]      : Número de threads a usar (padrão: 4)"
    echo "    [genes_list]   : Opcional: genes separados por vírgulas (e.g. ERG11,GSC1). Restringe a análise."
    exit 1
}

# Verifica se os argumentos obrigatórios foram fornecidos
[ -n "$GENES_FASTA" ] && [ -n "$RESULTS_DIR" ] && [ -n "$FASTQ_DIR" ] || print_usage

# Define o filtro de região para GATK/bcftools (apenas se GENES_LIST for fornecido)
REGION_FILTER=""
if [ -n "$GENES_LIST" ]; then
    echo "Aviso: A análise será RESTRITA aos genes: $GENES_LIST."
    # Transforma 'G1,G2,G3' no formato '-L G1 -L G2 -L G3' para GATK
    REGION_FILTER=$(echo "$GENES_LIST" | tr ',' '\n' | sed 's/^/-L /' | tr '\n' ' ')
fi


# --- FUNÇÕES DE LIMPEZA E REPORTING ---

# Função de limpeza de ficheiros temporários (definida aqui para ser usada no loop)
cleanup_temp_files() {
    local ISOLATE_DIR=$1
    local ISOLATE_NAME=$2
    # Variáveis temporárias (necessário redefini-las localmente)
    local TEMP_DECOMPOSED="$ISOLATE_DIR/$ISOLATE_NAME.decomposed.vcf"
    local TEMP_SNPS_VCF="$ISOLATE_DIR/$ISOLATE_NAME.snps.vcf.gz"
    local TEMP_INDELS_VCF="$ISOLATE_DIR/$ISOLATE_NAME.indels.vcf.gz"
    local TEMP_SNPS_FILTERED="$ISOLATE_DIR/$ISOLATE_NAME.snps.filtered.vcf.gz"
    local TEMP_INDELS_FILTERED="$ISOLATE_DIR/$ISOLATE_NAME.indels.filtered.vcf.gz"
    local PASS_POSITIONS_FILE="$ISOLATE_DIR/$ISOLATE_NAME.pass_positions.tmp"
    local UNPHASED_VCF="$ISOLATE_DIR/$ISOLATE_NAME${TEMP_VCF_SUFFIX}"

    # Limpeza dos ficheiros temporários intermédios
    rm -f "$TEMP_DECOMPOSED"
    rm -f "$TEMP_SNPS_VCF" "$TEMP_SNPS_VCF.tbi"
    rm -f "$TEMP_INDELS_VCF" "$TEMP_INDELS_VCF.tbi"
    rm -f "$TEMP_SNPS_FILTERED" "$TEMP_SNPS_FILTERED.tbi"
    rm -f "$TEMP_INDELS_FILTERED" "$TEMP_INDELS_FILTERED.tbi"
    rm -f "$PASS_POSITIONS_FILE" # Limpeza do ficheiro temporário
    rm -f ${UNPHASED_VCF}
    rm -f ${UNPHASED_VCF%.*}.tbi
}
# -------------------------------------


# ---- PHASE 0: Mapping FASTQs and generating BAMs ----

echo "### STEP 0: Mapping reads and generating sorted BAMs..."
echo "    -> Usando $THREADS threads."

if [ ! -f "${GENES_FASTA}.bwt" ]; then
    echo "    -> Indexando $GENES_FASTA para BWA..."
    bwa index "$GENES_FASTA"
fi

mkdir -p "$RESULTS_DIR"

for fq1 in "$FASTQ_DIR"/*_1*.fastq*; do
    [ -e "$fq1" ] || continue
    SAMPLE=$(basename "$fq1" | sed 's/_1.*//')
    fq2="$FASTQ_DIR/${SAMPLE}_2.fastq.gz"
    [ -f "$fq2" ] || fq2="$FASTQ_DIR/${SAMPLE}_2.fastq"
    if [ ! -f "$fq2" ]; then
        echo "    -> Aviso: Ficheiro R2 (par) não encontrado para $fq1. A saltar amostra $SAMPLE."
        continue
    fi
    smdir="$RESULTS_DIR/$SAMPLE"
    mkdir -p "$smdir"
    outbam="$smdir/${SAMPLE}.sorted.bam"
    if [ -s "$outbam" ]; then
        echo "    [$SAMPLE] BAM ordenado existe. A saltar alinhamento."
        continue
    fi
    echo "    [$SAMPLE] A alinhar e ordenar com BWA/SAMtools (Header @RG)..."
    # Adiciona o cabeçalho @RG adequado diretamente no comando bwa mem
    bwa mem -t "$THREADS" -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${RG_LB}\tPL:${RG_PL}\tPU:${RG_PU}" "$GENES_FASTA" "$fq1" "$fq2" |
        samtools view -bh -@ "$THREADS" - |
        samtools sort -@ "$THREADS" -o "$outbam"
    samtools index "$outbam"
    echo "    [$SAMPLE] BAM Concluído: $outbam"
done

# ---- STEP 1: Prepare reference indices for variant calling ----

echo "### STEP 1: Preparando índices de referência para GATK/BCFTOOLS..."

# 1.1. GATK requer um índice .dict para a referência FASTA.
if [ ! -f "${GENES_FASTA%.*}.dict" ]; then
    echo "    -> A criar índice .dict (Dicionário de Sequência) para GATK..."
    gatk CreateSequenceDictionary -R ${GENES_FASTA}
fi

# 1.2. Samtools/GATK também requer o índice .fai (FASTA Index)
if [ ! -f "${GENES_FASTA}.fai" ]; then
    echo "    -> A criar índice .fai (FASTA Index) para Samtools/pysam..."
    samtools faidx ${GENES_FASTA}
fi


# --- PREPARAÇÃO DO LOOP ---
ALL_ISOLATES=$(find "$RESULTS_DIR" -maxdepth 1 -mindepth 1 -type d -printf "%f\n" | grep -v '^\.$' | grep -v '^..$' | grep -v "^\.$")
NUM_ISOLATES=$(echo "$ALL_ISOLATES" | wc -l)

if [ "$NUM_ISOLATES" -eq 0 ]; then
    echo "Erro: Não foram encontrados subdiretórios de isolados (BAMs) em $RESULTS_DIR."
    exit 1
fi

echo "A processar $NUM_ISOLATES isolados encontrados em $RESULTS_DIR..."
CURRENT_COUNT=0

# ---- PHASE 2: Per-sample variant calling, phasing, and haplotype reconstruction (GATK -> WhatsHap -> Python) ----

echo "### STEP 2: Chamada, Filtragem e Faseamento (GATK, WhatsHap, Python)"
echo "------------------------------------------------------"

for ISOLATE_NAME in $ALL_ISOLATES; do
    CURRENT_COUNT=$((CURRENT_COUNT + 1))
    ISOLATE_DIR="$RESULTS_DIR/$ISOLATE_NAME"
    INPUT_BAM_TO_USE="$ISOLATE_DIR/$ISOLATE_NAME.sorted.bam" # BAM criado no STEP 0
    UNPHASED_VCF="$ISOLATE_DIR/$ISOLATE_NAME${TEMP_VCF_SUFFIX}"
    FILTERED_VCF="$ISOLATE_DIR/$ISOLATE_NAME${FILTERED_VCF_SUFFIX}" # VCF temporário filtrado
    PHASED_VCF="$ISOLATE_DIR/$ISOLATE_NAME${FINAL_VCF_SUFFIX}"
    IMBALANCE_REPORT="$ISOLATE_DIR/$ISOLATE_NAME${IMBALANCE_REPORT_SUFFIX}"
    ALL_VARIANTS_REPORT="$ISOLATE_DIR/$ISOLATE_NAME${ALL_VARIANTS_REPORT_SUFFIX}"

    # Arquivo temporário de posições PASS para reconciliação
    PASS_POSITIONS_FILE="$ISOLATE_DIR/$ISOLATE_NAME.pass_positions.tmp"

    # Variáveis temporárias para a filtragem:
    TEMP_DECOMPOSED="$ISOLATE_DIR/$ISOLATE_NAME.decomposed.vcf" # VCF (não comprimido)
    TEMP_SNPS_VCF="$ISOLATE_DIR/$ISOLATE_NAME.snps.vcf.gz"
    TEMP_INDELS_VCF="$ISOLATE_DIR/$ISOLATE_NAME.indels.vcf.gz"
    TEMP_SNPS_FILTERED="$ISOLATE_DIR/$ISOLATE_NAME.snps.filtered.vcf.gz"
    TEMP_INDELS_FILTERED="$ISOLATE_DIR/$ISOLATE_NAME.indels.filtered.vcf.gz"


    echo "[$CURRENT_COUNT/$NUM_ISOLATES] -> Processando isolado: $ISOLATE_NAME"

    # Verificação de ficheiro BAM de origem
    if [ ! -f "$INPUT_BAM_TO_USE" ]; then
        echo "    Aviso: Ficheiro BAM não encontrado em $INPUT_BAM_TO_USE. A saltar $ISOLATE_NAME."
        continue
    fi
    # ------------------------------------------------------------------------


    # --- 2.1. CHAMADA DE VARIANTES COM GATK (HaplotypeCaller) ---
    echo "    -> 2.1. A chamar variantes com GATK (cria $TEMP_VCF_SUFFIX)..."
    if [ -n "$REGION_FILTER" ]; then
        echo "    -> CRÍTICO: A usar filtro de região: $GENES_LIST"
    fi

    # Usamos o nome da pasta como --sample-name
    gatk HaplotypeCaller \
        -R ${GENES_FASTA} \
        -I ${INPUT_BAM_TO_USE} \
        -O ${UNPHASED_VCF} \
        ${REGION_FILTER} \
        --minimum-mapping-quality 30 \
        -ploidy 2 \
        --sample-name "${ISOLATE_NAME}" \
        --min-base-quality-score 20

    # Verifica se o GATK falhou antes de continuar
    if [ $? -ne 0 ]; then
        echo "    ERRO CRÍTICO GATK para o isolado $ISOLATE_NAME. A saltar faseamento/reconstrução."
        cleanup_temp_files "$ISOLATE_DIR" "$ISOLATE_NAME"
        continue
    fi

    # --- DEFINIÇÃO DOS LIMIARES DE FILTRAGEM (INCLUI QD MÍNIMO E NOVO FILTRO AF) ---
    SNPS_RPRS="-8.0"
    SNPS_BQRS="-8.0"
    SNPS_AD_ALT="10"
    SNPS_DP_TOTAL="20"
    SNPS_QD="2.0"

    INDELS_RPRS="-20.0"
    INDELS_BQRS="-8.0"
    INDELS_AD_ALT="10"
    INDELS_DP_TOTAL="20"
    INDELS_QD="10.0" # Alterado de 25.0 para 10.0

    # NOVO: Limites de Profundidade Alélica (AF) para Heterozigotos (AF deve ser entre 0.4 e 0.6)
    AF_MIN_BALANCE="0.4"
    AF_MAX_BALANCE="0.6"

    # VALOR PADRÃO para Profundidade Máxima (será substituído pelo cálculo dinâmico)
    MAX_DP_FILTER="1000"

    # --- CÁLCULO DINÂMICO DO MAX_DP_FILTER (ROBUSTO) ---
    echo "    -> Calculando dinamicamente o MAX_DP_FILTER com base na DP média do isolado..."

    DP_VALUES=$(bcftools query -f '[%DP]\n' ${UNPHASED_VCF} -i 'FORMAT/DP>0' 2>/dev/null)
    NUM_VARIANTS=$(echo "$DP_VALUES" | wc -l)

    if [ "$NUM_VARIANTS" -gt 0 ]; then
        # Adicionado || true ao bc para evitar falha se o input for inválido/vazio
        DP_MEAN_FLOAT=$(echo "$DP_VALUES" | awk 'BEGIN {sum=0} {sum+=$1} END {printf "%.3f", sum/NR}' || true)

        if [ -n "$DP_MEAN_FLOAT" ]; then
            NEW_MAX_DP_FLOAT=$(echo "scale=3; $DP_MEAN_FLOAT * $DP_MULTIPLIER" | bc)
            NEW_MAX_DP_FILTER=$(printf "%.0f\n" $NEW_MAX_DP_FLOAT)

            MIN_LIMIT=$((SNPS_DP_TOTAL * 2))

            if [ "$NEW_MAX_DP_FILTER" -gt "$MIN_LIMIT" ] 2>/dev/null; then
                MAX_DP_FILTER="$NEW_MAX_DP_FILTER"
                echo "        -> DP Média calculada: $(printf "%.2f" "$DP_MEAN_FLOAT")X"
                echo "        -> Limite Máximo DP (Arredondado): ${MAX_DP_FILTER}X (Dinâmico)"
            else
                echo "        -> Aviso: Cálculo resultou num limite muito baixo. Mantendo valor padrão: ${MAX_DP_FILTER}X."
            fi
        else
            echo "        -> Aviso: Falha no cálculo DP Média. Mantendo valor padrão: ${MAX_DP_FILTER}X."
        fi
    else
        echo "        -> Aviso: Nenhuma variante válida (DP>0) para calcular DP Média. Mantendo valor padrão: ${MAX_DP_FILTER}X."
    fi
    # -------------------------------------------------------------------------


    # --- DEFINIÇÃO DAS EXPRESSÕES DE FILTRAGEM ---

    # 0. NON_RANKSUM_FILTER (Condições de Profundidade e Qualidade Mínima)
    NON_RANKSUM_SNPS="FMT/AD[0:1] > ${SNPS_AD_ALT} & FMT/DP[0] > ${SNPS_DP_TOTAL} & FMT/DP[0] < ${MAX_DP_FILTER} & INFO/QD > ${SNPS_QD}"
    NON_RANKSUM_INDELS="FMT/AD[0:1] > ${INDELS_AD_ALT} & FMT/DP[0] > ${INDELS_DP_TOTAL} & FMT/DP[0] < ${MAX_DP_FILTER} & INFO/QD > ${INDELS_QD}"

    # 1. PATH A: PERFECT 1/1 (Bypass RankSum - Condição: GT=1/1 E AD_REF=0)
    PERFECT_1_1_SNPS="FMT/GT=\"1/1\" & FMT/AD[0:0]=0 & (${NON_RANKSUM_SNPS})"
    PERFECT_1_1_INDELS="FMT/GT=\"1/1\" & FMT/AD[0:0]=0 & (${NON_RANKSUM_INDELS})"

    # 2. PATH B: STANDARD QUALITY CHECK (Aplica todos os RankSums + AF Balance)
    STRICT_Q_SNPS="(ReadPosRankSum > ${SNPS_RPRS} & BaseQRankSum > ${SNPS_BQRS})"
    STRICT_Q_INDELS="(ReadPosRankSum > ${INDELS_RPRS} & BaseQRankSum > ${INDELS_BQRS})"

    STANDARD_AF_GT_CHECK="(FMT/GT=\"1/1\" | (FMT/GT=\"0/1\" & (FMT/AD[0:1] / FMT/DP[0]) >= ${AF_MIN_BALANCE} & (FMT/AD[0:1] / FMT/DP[0]) <= ${AF_MAX_BALANCE}))"

    STANDARD_PASS_SNPS="(${NON_RANKSUM_SNPS} & ${STRICT_Q_SNPS} & ${STANDARD_AF_GT_CHECK})"
    STANDARD_PASS_INDELS="(${NON_RANKSUM_INDELS} & ${STRICT_Q_INDELS} & ${STANDARD_AF_GT_CHECK})"

    # 3. FILTRO FINAL: PATH A (PERFECT 1/1) OU PATH B (STANDARD PASS)
    SNPS_FINAL_FILTER="(${PERFECT_1_1_SNPS}) | (${STANDARD_PASS_SNPS})"
    INDELS_FINAL_FILTER="(${PERFECT_1_1_INDELS}) | (${STANDARD_PASS_INDELS})"

    # 4. Filtro para o RELATÓRIO DE IMBALANCE (0/1 Desequilibrado que passe a Qualidade Básica)
    IMBALANCE_AF_CHECK_FILTER="FMT/GT=\"0/1\" & ((FMT/AD[0:1] / FMT/DP[0]) < ${AF_MIN_BALANCE} | (FMT/AD[0:1] / FMT/DP[0]) > ${AF_MAX_BALANCE})"
    # -------------------------------------------------------------------------


    # --- 2.1.5. FILTRAGEM POR SUPORTE DE READS (AD) E QUALIDADE (RankSums + QD) ---
    echo "    -> 2.1.5. A aplicar filtragem robusta (DP < ${MAX_DP_FILTER}X Dinâmico)..."
    echo "    -> Novo filtro de caminho duplo implementado para aceitar 1/1 perfeitos (AD_REF=0)."

    COUNT_BEFORE=$(bcftools view -H ${UNPHASED_VCF} 2>/dev/null | wc -l)
    echo "        -> Variantes antes da filtragem: ${COUNT_BEFORE} (linhas de dados)"

    if [ "$BYPASS_FILTER" = "true" ]; then
        echo "        [DEBUG MODO] FILTRAGEM IGNORADA (BYPASS). Apenas copiando..."
        # Copia e garante BGZF e indexação
        bcftools view -O z -o "$FILTERED_VCF" ${UNPHASED_VCF}
    else
        # 2.1.5.a. Decompor blocos multialélicos e normalizar
        echo "        -> 2.1.5.a. Decompondo/Normalizando variantes (bcftools norm)..."
        bcftools norm -f ${GENES_FASTA} -m -any -O v -o "$TEMP_DECOMPOSED" ${UNPHASED_VCF}

        if [ $? -ne 0 ]; then
            echo "    ERRO CRÍTICO BCFTOOLS NORM: Falha na decomposição para o isolado $ISOLATE_NAME. A saltar."
            cleanup_temp_files "$ISOLATE_DIR" "$ISOLATE_NAME"
            continue
        fi

        # -------------------------------------------------------------------------
        # --- 2.1.2. GERAÇÃO DO RELATÓRIO COMPLETO DE FILTRAGEM (PRIMEIRA PASSAGEM - FAIL REASONS) ---
        echo "    -> 2.1.2. A gerar Relatório Completo de Filtragem ($ALL_VARIANTS_REPORT_SUFFIX) com motivos de falha (Status PENDING)..."

        # Expressões de avaliação do filtro
        EVAL_SNPS_Q_FILTER="${NON_RANKSUM_SNPS}"
        EVAL_INDELS_Q_FILTER="${NON_RANKSUM_INDELS}"
        EVAL_AF_BALANCE_FILTER="FMT/GT=\"0/1\" & (FMT/AD[0:1] / FMT/DP[0]) >= ${AF_MIN_BALANCE} & (FMT/AD[0:1] / FMT/DP[0]) <= ${AF_MAX_BALANCE}"

        # Cabeçalho do Relatório
        echo -e "# Relatório Completo de Avaliação de Variantes\n# ISOLADO: $ISOLATE_NAME" > "$ALL_VARIANTS_REPORT"
        echo -e "# Definições de Limiares:\n# Min DP Total: $SNPS_DP_TOTAL, Max DP Total: $MAX_DP_FILTER, Min AD Alt: $SNPS_AD_ALT, QD Min SNP/Indel: $SNPS_QD/$INDELS_QD, RankSum Min SNP/Indel: $SNPS_RPRS/$INDELS_RPRS, AF Balance: $AF_MIN_BALANCE-$AF_MAX_BALANCE" >> "$ALL_VARIANTS_REPORT"
        echo -e "# CHROM\tPOS\tREF\tALT\tTIPO\tQUAL\tDP\tAD_REF,AD_ALT\tGT\tAF_ALT_CALCULADO\tSTATUS_FINAL\tFAIL_REASONS" >> "$ALL_VARIANTS_REPORT"

        # Função para avaliar e reportar uma variante (mantida a lógica do script anterior)
        report_variant() {
            local TYPE=$1
            local Q_FILTER=$2

            # Filtra variantes GT=1/1 ou GT=0/1 para serem relevantes, ignorando 0/0.
            # Adicionado || true ao final do pipe para evitar falha do 'set -e' se não houver output
            bcftools view -H -i 'FMT/GT="1/1" | FMT/GT="0/1"' -v $TYPE "$TEMP_DECOMPOSED" 2>/dev/null | \
            while read -r LINE; do
                # 1. Extrair os dados da variante
                CHROM=$(echo "$LINE" | awk '{print $1}')
                POS=$(echo "$LINE" | awk '{print $2}')
                REF=$(echo "$LINE" | awk '{print $4}')
                ALT=$(echo "$LINE" | awk '{print $5}')
                QUAL=$(echo "$LINE" | awk '{print $6}')

                FORMAT_FIELD=$(echo "$LINE" | awk '{print $9}')
                SAMPLE_FIELD=$(echo "$LINE" | awk '{print $10}')

                GT_INDEX=$(echo "$FORMAT_FIELD" | tr ':' '\n' | grep -n '^GT$' | cut -d: -f1)
                AD_INDEX=$(echo "$FORMAT_FIELD" | tr ':' '\n' | grep -n '^AD$' | cut -d: -f1)
                DP_INDEX=$(echo "$FORMAT_FIELD" | tr ':' '\n' | grep -n '^DP$' | cut -d: -f1)

                if [ -z "$GT_INDEX" ] || [ -z "$AD_INDEX" ] || [ -z "$DP_INDEX" ]; then
                    echo -e "$CHROM\t$POS\t$REF\t$ALT\t$TYPE\t$QUAL\tN/A\tN/A\tN/A\tN/A\tPENDING\tMISSING_FORMAT_FIELD" >> "$ALL_VARIANTS_REPORT"
                    continue
                fi

                GT=$(echo "$SAMPLE_FIELD" | cut -d: -f$GT_INDEX)
                AD=$(echo "$SAMPLE_FIELD" | cut -d: -f$AD_INDEX)
                DP=$(echo "$SAMPLE_FIELD" | cut -d: -f$DP_INDEX)
                AD_ALT=$(echo "$AD" | cut -d, -f2)

                AF_CALCULADO="N/A"
                # Adicionado 2>/dev/null ao bc para suprimir erros se DP for zero/ponto
                if [ -n "$DP" ] && [ "$DP" -ne 0 ] && [ "$AD_ALT" != "." ]; then
                    AF_CALCULADO=$(echo "scale=3; $AD_ALT / $DP" | bc 2>/dev/null)
                fi

                # 2. Avaliar os critérios de FALHA (para gerar as RAZÕES)
                # Adicionado || true para garantir que o bcftools filter não cause a saída do script
                IS_Q_PASS=$(bcftools filter -i "${Q_FILTER}" -r $CHROM:$POS-$POS "$TEMP_DECOMPOSED" 2>/dev/null | grep -v '^#' | wc -l || true)

                IS_AF_BALANCE_PASS=1
                if [ "$GT" = "0/1" ]; then
                    # Adicionado || true para garantir que o bcftools filter não cause a saída do script
                    IS_AF_BALANCE_PASS=$(bcftools filter -i "${EVAL_AF_BALANCE_FILTER}" -r $CHROM:$POS-$POS "$TEMP_DECOMPOSED" 2>/dev/null | grep -v '^#' | wc -l || true)
                fi

                # 3. Determinar RAZÕES DE FALHA
                FINAL_STATUS="PENDING"
                FAIL_REASONS=""

                IS_FAILED_Q=0
                if [ "$IS_Q_PASS" -eq 0 ]; then
                    FAIL_REASONS="${FAIL_REASONS}Q_FAIL(DP_AD_QD_MaxDP);"
                    IS_FAILED_Q=1
                fi

                IS_FAILED_AF=0
                if [ "$GT" = "0/1" ] && [ "$IS_AF_BALANCE_PASS" -eq 0 ]; then
                    FAIL_REASONS="${FAIL_REASONS}AF_IMBALANCE(0.4-0.6);"
                    IS_FAILED_AF=1
                fi

                if [ "$IS_FAILED_Q" -eq 1 ] || [ "$IS_FAILED_AF" -eq 1 ]; then
                    FINAL_STATUS="FAIL_PRELIMINAR"
                fi

                # 4. Escrever a linha no relatório
                echo -e "$CHROM\t$POS\t$REF\t$ALT\t$TYPE\t$QUAL\t$DP\t$AD\t$GT\t$AF_CALCULADO\t$FINAL_STATUS\t$FAIL_REASONS" >> "$ALL_VARIANTS_REPORT"

            done
        }

        # Gerar relatório para SNPs e depois para Indels
        # O pipe da função report_variant é seguro (com || true)
        report_variant "snps" "${EVAL_SNPS_Q_FILTER}"
        report_variant "indels" "${EVAL_INDELS_Q_FILTER}"

        REPORT_COUNT_TOTAL=$(wc -l < "$ALL_VARIANTS_REPORT" 2>/dev/null || echo 0)

        if [ "$REPORT_COUNT_TOTAL" -ge 4 ]; then
            REPORT_DATA_LINES=$((REPORT_COUNT_TOTAL - 3))
            echo "    -> ${ALL_VARIANTS_REPORT_SUFFIX} gerado com $REPORT_DATA_LINES variantes avaliadas."
        else
            echo "    -> Aviso: ${ALL_VARIANTS_REPORT_SUFFIX} gerado, mas encontrou poucas variantes para reportar."
        fi
        # -------------------------------------------------------------------------

        # 2.1.5.b. Filtragem de SNPs (Path A OR Path B)
        echo "        -> 2.1.5.b. Selecionando e filtrando SNPs..."
        bcftools view -v snps -O z -o "$TEMP_SNPS_VCF" "$TEMP_DECOMPOSED"
        tabix -p vcf "$TEMP_SNPS_VCF"

        bcftools filter -i "${SNPS_FINAL_FILTER}" -o "$TEMP_SNPS_FILTERED" -O z "$TEMP_SNPS_VCF"
        tabix -p vcf "$TEMP_SNPS_FILTERED"

        # 2.1.5.c. Filtragem de Indels (Path A OR Path B)
        echo "        -> 2.1.5.c. Selecionando e filtrando Indels..."
        bcftools view -v indels -O z -o "$TEMP_INDELS_VCF" "$TEMP_DECOMPOSED"
        tabix -p vcf "$TEMP_INDELS_VCF"

        bcftools filter -i "${INDELS_FINAL_FILTER}" -o "$TEMP_INDELS_FILTERED" -O z "$TEMP_INDELS_VCF"
        tabix -p vcf "$TEMP_INDELS_FILTERED"

        # 2.1.5.d. Recombinação (cria o VCF final PASS)
        echo "        -> 2.1.5.d. Recombinando SNPs e Indels filtrados. Output final: ${FILTERED_VCF}..."

        VCFS_TO_CONCAT=""
        if [ -s "$TEMP_SNPS_FILTERED" ]; then VCFS_TO_CONCAT="$VCFS_TO_CONCAT $TEMP_SNPS_FILTERED"; fi
        if [ -s "$TEMP_INDELS_FILTERED" ]; then VCFS_TO_CONCAT="$VCFS_TO_CONCAT $TEMP_INDELS_FILTERED"; fi

        if [ -z "$VCFS_TO_CONCAT" ]; then
            echo "        -> Aviso: Nenhuma variante restante após a filtragem. A criar VCF vazio válido."
            # Cria um VCF BGZF válido, mas vazio, mantendo o cabeçalho
            bcftools view -h ${UNPHASED_VCF} 2>/dev/null | bgzip -c > "$FILTERED_VCF"
        else
            bcftools concat $VCFS_TO_CONCAT -a -O z -o "$FILTERED_VCF"
        fi

        # --- 2.1.7. CRÍTICO: RECONCILIAÇÃO DO RELATÓRIO TXT COM O VCF FILTRADO ---
        echo "    -> 2.1.7. RECONCILIAÇÃO: A corrigir o status PASS/FAIL no relatório TXT..."

        # 1. Extrai as posições (CHROM e POS) de todas as variantes que *realmente* passaram
        bcftools query -f '%CHROM\t%POS\n' ${FILTERED_VCF} 2>/dev/null > "$PASS_POSITIONS_FILE"

        # 2. Reconciliação com AWK (mantida a lógica exata para a reconciliação)
        TEMP_REPORT="$ALL_VARIANTS_REPORT.tmp"

        awk -v OFS='\t' '
            NR==FNR {
                if (FNR > 0) { pass_positions[$1 ":" $2] = 1; }
                next;
            }
            FNR==1 { print; next; } # Imprime o cabeçalho do relatório

            {
                key = $1 ":" $2;

                if (key in pass_positions) {
                    $11 = "PASS";  # Coluna de STATUS_FINAL
                    $12 = "";      # Coluna de FAIL_REASONS (limpa)
                }
                else if ($11 == "FAIL_PRELIMINAR") {
                    $11 = "FAIL";
                }
                else if ($11 == "PENDING") {
                    $11 = "FAIL";
                    $12 = "FILTER_FAIL_OTHER";
                }

                print;
            }
        ' "$PASS_POSITIONS_FILE" "$ALL_VARIANTS_REPORT" > "$TEMP_REPORT"

        mv "$TEMP_REPORT" "$ALL_VARIANTS_REPORT"
        echo "    -> Status do relatório corrigido: $ALL_VARIANTS_REPORT agora corresponde ao VCF final."
        # -------------------------------------------------------------------------


        # 2.1.5.e. CRÍTICO: GERAÇÃO DO RELATÓRIO DE VARIANTE DESEQUILIBRADAS
        echo "        -> 2.1.5.e. A gerar Relatório de Variantes Desequilibradas (0/1 imbalanced)..."

        # EXPRESSÃO COMPLETA PARA RELATÓRIO DE IMBALANCE
        SNPS_IMBALANCE_REPORT_FILTER="(${NON_RANKSUM_SNPS}) & (${IMBALANCE_AF_CHECK_FILTER})"
        INDELS_IMBALANCE_REPORT_FILTER="(${NON_RANKSUM_INDELS}) & (${IMBALANCE_AF_CHECK_FILTER})"

        # Cria o cabeçalho do relatório
        echo -e "# Relatório de Variantes Desequilibradas de Alta Qualidade (FMT/GT=0/1, mas AF fora de 0.4-0.6)\n# CHROM\tPOS\tREF\tALT\tQUAL\tAD\tDP\tGT\tAF_CALCULADO" > "$IMBALANCE_REPORT"

        # Query de SNPs desequilibrados (usa o VCF decomposed/normalized)
        # Adicionado || true ao final para evitar falha do 'set -e' se não houver output
        bcftools filter -i "${SNPS_IMBALANCE_REPORT_FILTER}" "$TEMP_DECOMPOSED" 2>/dev/null | \
        bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%AD\t%DP\t%GT\t%AD[0:1]/%DP[0]\n]' >> "$IMBALANCE_REPORT" || true

        # Query de INDELs desequilibrados
        # Adicionado || true ao final para evitar falha do 'set -e' se não houver output
        bcftools filter -i "${INDELS_IMBALANCE_REPORT_FILTER}" "$TEMP_DECOMPOSED" 2>/dev/null | \
        bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%AD\t%DP\t%GT\t%AD[0:1]/%DP[0]\n]' >> "$IMBALANCE_REPORT" || true

        REPORT_COUNT=$(wc -l < "$IMBALANCE_REPORT" 2>/dev/null || echo 1)
        echo "        -> ${IMBALANCE_REPORT_SUFFIX} gerado com $((REPORT_COUNT - 1)) variante(s) suspeita(s)."

    fi

    # --- 2.1.8. DEBUG: Listar as primeiras 5 variantes de cada gene/contig ---
    echo "    -> 2.1.8. DEBUG: A listar as primeiras 5 variantes detectadas para CADA contig/gene..."

    # Indexar o VCF filtrado final, pois o WhatsHap precisa.
    echo "    -> A indexar o VCF filtrado para WhatsHap..."
    tabix -p vcf ${FILTERED_VCF}

    # Adicionado || true ao find para evitar falha se o VCF estiver vazio (e bcftools não encontrar contigs)
    ALL_CONTIGS=$(bcftools view -H ${FILTERED_VCF} 2>/dev/null | awk '{print $1}' | sort -u || true)

    if [ -z "$ALL_CONTIGS" ]; then
        echo "        -> Aviso: Nenhuma variante encontrada no VCF filtrado. (Verifique os limiares)."
    else
        echo "        --- INÍCIO DO RELATÓRIO DE VARIANTE POR CONTIG (MAX 20 LINHAS CADA) ---"
        for CONTIG in $ALL_CONTIGS; do
            echo ""
            echo "        == CONTIG/GENE: $CONTIG (Primeiras 20 Variantes) =="
            bcftools query -r $CONTIG -f '[%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%AD\t%DP\t%GT\t%AD[0:1]/%DP[0]\n]' ${FILTERED_VCF} | \
            awk 'BEGIN {OFS="\t"} {$6="PASS"; print}' | head -n 20
        done
        echo "        ---------------------------------------------------"
    fi
    # ----------------------------------------------------------------------


    # --- 2.2. FASEAMENTO (PHASING) COM WHATSHAP ---
    echo "    -> 2.2. A fasear variantes com WhatsHap (cria $FINAL_VCF_SUFFIX)..."

    whatshap polyphase \
        --ploidy 2 \
        --reference ${GENES_FASTA} \
        --distrust-genotypes \
        ${FILTERED_VCF} \
        ${INPUT_BAM_TO_USE} \
        -o ${PHASED_VCF}

    # Verifica se o WhatsHap falhou
    if [ $? -ne 0 ]; then
        echo "    ERRO CRÍTICO WHATSHAP: Falha no faseamento para o isolado $ISOLATE_NAME. A saltar reconstrução."
        cleanup_temp_files "$ISOLATE_DIR" "$ISOLATE_NAME"
        continue
    fi

    # Limpeza dos ficheiros VCF filtrados (que foram input para o WhatsHap)
    rm -f ${FILTERED_VCF}
    rm -f ${FILTERED_VCF%.*}.tbi

    # --- 2.2.6. INDEXAÇÃO DO VCF FASEADO (Para o script Python) ---
    echo "    -> 2.2.6. A indexar o VCF faseado com tabix..."
    tabix -p vcf ${PHASED_VCF}

    # --- 2.3. RECONSTRUÇÃO DE HAPLÓTIPOS COM PYTHON ---
    echo "    -> 2.3. A reconstruir haplótipos (Python)..."
    python ${HAPLOTYPE_SCRIPT} --reconstruct_single ${GENES_FASTA} ${RESULTS_DIR} ${ISOLATE_NAME}

    # --- 2.4. LIMPEZA FINAL ---
    cleanup_temp_files "$ISOLATE_DIR" "$ISOLATE_NAME"

    echo "[$CURRENT_COUNT/$NUM_ISOLATES] -> $ISOLATE_NAME CONCLUÍDO."
    echo "" # Linha vazia para clareza
done

# ---- PHASE 3: Consolidação final ----

echo "### STEP 3: Consolidação final"
echo "------------------------------------------------------"
echo "3. A consolidar todos os haplótipos num FASTA por gene..."
echo "------------------------------------------------------"

# Adicionado || true caso o script python falhe por falta de output (mas o output devia ser gerado)
python ${HAPLOTYPE_SCRIPT} --consolidate_all ${RESULTS_DIR} || true

echo "------------------------------------------------------"
echo "PIPELINE CONCLUÍDA! Os ficheiros FASTA finais por gene estão em $RESULTS_DIR."
echo "------------------------------------------------------"
