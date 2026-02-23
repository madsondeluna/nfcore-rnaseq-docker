# Pipeline de RNA-Seq com nf-core/rnaseq e Docker

Guia completo para executar analise de RNA-Seq -- desde a instalacao do Docker ate a extracao de genes diferencialmente expressos (DEGs) -- utilizando o pipeline [nf-core/rnaseq](https://nf-co.re/rnaseq) com containers Docker.

---

## Sumario

1. [Visao Geral](#visao-geral)
2. [Pre-requisitos](#pre-requisitos)
3. [Etapa 1 -- Instalacao e Configuracao do Docker](#etapa-1----instalacao-e-configuracao-do-docker)
4. [Etapa 2 -- Instalacao do Nextflow](#etapa-2----instalacao-do-nextflow)
5. [Etapa 3 -- Organizacao dos Dados](#etapa-3----organizacao-dos-dados)
6. [Etapa 4 -- Preparacao do Samplesheet](#etapa-4----preparacao-do-samplesheet)
7. [Etapa 5 -- Obtencao do Genoma de Referencia](#etapa-5----obtencao-do-genoma-de-referencia)
8. [Etapa 6 -- Execucao do Pipeline nf-core/rnaseq](#etapa-6----execucao-do-pipeline-nf-corernaseq)
9. [Etapa 7 -- Interpretacao dos Resultados](#etapa-7----interpretacao-dos-resultados)
10. [Etapa 8 -- Analise de Expressao Diferencial (DEGs)](#etapa-8----analise-de-expressao-diferencial-degs)
11. [Resolucao de Problemas](#resolucao-de-problemas)
12. [Referencias](#referencias)

---

## Visao Geral

Este protocolo utiliza o pipeline **nf-core/rnaseq** (versao 3.18.0 ou superior) para processar dados brutos de RNA-Seq (arquivos FASTQ) e gerar matrizes de contagem de genes. Posteriormente, utiliza o pipeline **nf-core/differentialabundance** para identificar genes diferencialmente expressos (DEGs).

### O que o pipeline nf-core/rnaseq faz

O pipeline executa as seguintes etapas automaticamente:

1. **Controle de qualidade dos reads brutos** -- FastQC
2. **Trimagem de adaptadores** -- Trim Galore / fastp
3. **Alinhamento ao genoma de referencia** -- STAR ou HISAT2
4. **Quantificacao de transcritos** -- Salmon, RSEM ou featureCounts
5. **Controle de qualidade pos-alinhamento** -- RSeQC, Qualimap, dupRadar
6. **Relatorio consolidado** -- MultiQC

### Fluxo geral

```
FASTQ brutos
    |
    v
[Controle de Qualidade - FastQC]
    |
    v
[Trimagem - Trim Galore]
    |
    v
[Alinhamento - STAR/HISAT2]
    |
    v
[Quantificacao - Salmon/featureCounts]
    |
    v
[Relatorio - MultiQC]
    |
    v
Matrizes de contagem
    |
    v
[Analise Diferencial - DESeq2/edgeR]
    |
    v
Lista de DEGs
```

---

## Pre-requisitos

| Requisito | Minimo | Recomendado |
|-----------|--------|-------------|
| RAM | 8 GB | 16 GB ou mais |
| Disco | 50 GB livres | 200 GB ou mais |
| CPU | 4 nucleos | 8 nucleos ou mais |
| Sistema Operacional | Linux, macOS ou Windows (via WSL2) | Linux |
| Java | versao 11 ou superior | versao 17 |
| Conexao com internet | necessaria para download de containers e genomas | -- |

---

## Etapa 1 -- Instalacao e Configuracao do Docker

O Docker permite executar softwares em containers isolados, garantindo que todas as dependencias estejam corretas e que o ambiente seja reprodutivel.

### 1.1 Instalacao no Linux (Ubuntu/Debian)

```bash
# Atualizar os pacotes do sistema
sudo apt-get update

# Instalar dependencias necessarias
sudo apt-get install -y \
    ca-certificates \
    curl \
    gnupg \
    lsb-release

# Adicionar a chave GPG oficial do Docker
sudo mkdir -p /etc/apt/keyrings
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | \
    sudo gpg --dearmor -o /etc/apt/keyrings/docker.gpg

# Adicionar o repositorio do Docker
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.gpg] \
  https://download.docker.com/linux/ubuntu \
  $(lsb_release -cs) stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null

# Instalar o Docker Engine
sudo apt-get update
sudo apt-get install -y docker-ce docker-ce-cli containerd.io docker-compose-plugin
```

### 1.2 Instalacao no macOS

1. Acessar o site oficial: https://www.docker.com/products/docker-desktop/
2. Baixar o **Docker Desktop** para macOS (Apple Silicon ou Intel, conforme seu hardware).
3. Abrir o arquivo `.dmg` baixado.
4. Arrastar o icone do Docker para a pasta **Applications**.
5. Abrir o Docker Desktop a partir do Launchpad ou da pasta Applications.
6. Aguardar a inicializacao completa (o icone da baleia na barra de menu ficara estavel).

### 1.3 Instalacao no Windows

1. Acessar https://www.docker.com/products/docker-desktop/
2. Baixar o **Docker Desktop** para Windows.
3. Executar o instalador e seguir as instrucoes na tela.
4. **Importante:** habilitar o WSL2 quando solicitado durante a instalacao.
5. Reiniciar o computador se solicitado.
6. Abrir o Docker Desktop e aguardar a inicializacao.

### 1.4 Configuracao pos-instalacao (Linux)

Por padrao no Linux, o Docker exige permissoes de superusuario. Para executar sem `sudo`:

```bash
# Adicionar seu usuario ao grupo docker
sudo usermod -aG docker $USER

# Aplicar a mudanca (necessario logout/login, ou executar)
newgrp docker
```

### 1.5 Verificar se o Docker esta funcionando

```bash
# Verificar a versao instalada
docker --version

# Executar o container de teste
docker run hello-world
```

Se a mensagem "Hello from Docker!" aparecer, a instalacao foi bem-sucedida.

### 1.6 Configurar recursos do Docker

Para analises de RNA-Seq, e recomendavel alocar recursos adequados ao Docker:

- **Docker Desktop (macOS/Windows):** abrir Docker Desktop > Settings > Resources > ajustar CPU, memoria e disco conforme necessidade.
- **Linux:** o Docker utiliza os recursos do host diretamente, sem necessidade de configuracao adicional.

---

## Etapa 2 -- Instalacao do Nextflow

O Nextflow e o motor de workflow que executa o pipeline nf-core/rnaseq. Ele orquestra os containers Docker automaticamente.

### 2.1 Verificar se o Java esta instalado

```bash
java -version
```

Se o Java nao estiver instalado ou for inferior a versao 11:

```bash
# Ubuntu/Debian
sudo apt-get install -y default-jdk

# macOS (com Homebrew)
brew install openjdk@17
```

### 2.2 Instalar o Nextflow

```bash
# Baixar e instalar o Nextflow
curl -s https://get.nextflow.io | bash

# Tornar o binario executavel
chmod +x nextflow

# Mover para um diretorio no PATH do sistema
sudo mv nextflow /usr/local/bin/
```

### 2.3 Verificar a instalacao

```bash
nextflow -version
```

### 2.4 (Opcional) Instalar as ferramentas de linha de comando do nf-core

As ferramentas nf-core nao sao obrigatorias para rodar o pipeline, mas facilitam tarefas como baixar pipelines e gerar samplesheets.

```bash
pip install nf-core
```

---

## Etapa 3 -- Organizacao dos Dados

### 3.1 Estrutura de diretorios recomendada

Criar a seguinte estrutura antes de iniciar a analise:

```
projeto_rnaseq/
├── data/
│   ├── raw_fastq/          # Arquivos FASTQ brutos aqui
│   │   ├── amostra1_R1.fastq.gz
│   │   ├── amostra1_R2.fastq.gz
│   │   ├── amostra2_R1.fastq.gz
│   │   ├── amostra2_R2.fastq.gz
│   │   └── ...
│   └── reference/           # Genoma de referencia
│       ├── genome.fa
│       └── genes.gtf
├── samplesheet.csv           # Samplesheet de entrada
├── results/                  # Resultados do pipeline
└── work/                     # Diretorio de trabalho do Nextflow
```

### 3.2 Criar os diretorios

```bash
mkdir -p projeto_rnaseq/{data/{raw_fastq,reference},results}
cd projeto_rnaseq
```

### 3.3 Transferir os arquivos FASTQ

Copiar os arquivos FASTQ brutos para o diretorio `data/raw_fastq/`. Os arquivos devem estar compactados com gzip (extensao `.fastq.gz` ou `.fq.gz`).

```bash
# Exemplo: copiar de um HD externo
cp /media/hd_externo/fastq/*.fastq.gz data/raw_fastq/

# Exemplo: baixar de um servidor via scp
scp usuario@servidor:/caminho/dos/dados/*.fastq.gz data/raw_fastq/
```

---

## Etapa 4 -- Preparacao do Samplesheet

O samplesheet e um arquivo CSV que informa ao pipeline quais amostras processar, onde estao os arquivos FASTQ e qual a orientacao da biblioteca (strandedness).

### 4.1 Formato do samplesheet

O arquivo deve conter exatamente 4 colunas, separadas por virgula, com cabecalho:

| Coluna | Descricao |
|--------|-----------|
| `sample` | Nome da amostra (identificador unico por amostra biologica) |
| `fastq_1` | Caminho completo para o arquivo FASTQ R1 (forward) |
| `fastq_2` | Caminho completo para o arquivo FASTQ R2 (reverse). Deixar vazio para single-end |
| `strandedness` | Orientacao da biblioteca: `unstranded`, `forward`, `reverse` ou `auto` |

### 4.2 Exemplo de samplesheet para dados paired-end

Criar o arquivo `samplesheet.csv`:

```csv
sample,fastq_1,fastq_2,strandedness
CONTROLE_REP1,/caminho/completo/data/raw_fastq/CTRL_1_R1.fastq.gz,/caminho/completo/data/raw_fastq/CTRL_1_R2.fastq.gz,auto
CONTROLE_REP2,/caminho/completo/data/raw_fastq/CTRL_2_R1.fastq.gz,/caminho/completo/data/raw_fastq/CTRL_2_R2.fastq.gz,auto
CONTROLE_REP3,/caminho/completo/data/raw_fastq/CTRL_3_R1.fastq.gz,/caminho/completo/data/raw_fastq/CTRL_3_R2.fastq.gz,auto
TRATAMENTO_REP1,/caminho/completo/data/raw_fastq/TREAT_1_R1.fastq.gz,/caminho/completo/data/raw_fastq/TREAT_1_R2.fastq.gz,auto
TRATAMENTO_REP2,/caminho/completo/data/raw_fastq/TREAT_2_R1.fastq.gz,/caminho/completo/data/raw_fastq/TREAT_2_R2.fastq.gz,auto
TRATAMENTO_REP3,/caminho/completo/data/raw_fastq/TREAT_3_R1.fastq.gz,/caminho/completo/data/raw_fastq/TREAT_3_R2.fastq.gz,auto
```

### 4.3 Exemplo para dados single-end

```csv
sample,fastq_1,fastq_2,strandedness
CONTROLE_REP1,/caminho/completo/data/raw_fastq/CTRL_1.fastq.gz,,auto
CONTROLE_REP2,/caminho/completo/data/raw_fastq/CTRL_2.fastq.gz,,auto
TRATAMENTO_REP1,/caminho/completo/data/raw_fastq/TREAT_1.fastq.gz,,auto
TRATAMENTO_REP2,/caminho/completo/data/raw_fastq/TREAT_2.fastq.gz,,auto
```

### 4.4 Sobre a coluna strandedness

| Valor | Quando usar |
|-------|-------------|
| `auto` | Na duvida, usar este. O pipeline detecta automaticamente usando Salmon |
| `unstranded` | Quando a biblioteca nao preserva informacao de fita (ex.: TruSeq unstranded) |
| `forward` | Biblioteca strand-specific no sentido forward (ex.: Ligation, dUTP - SMARTer) |
| `reverse` | Biblioteca strand-specific no sentido reverse (ex.: TruSeq Stranded, maioria dos kits Illumina atuais) |

### 4.5 Regras importantes

- Usar **caminhos absolutos** (completos) para os arquivos FASTQ.
- Se uma mesma amostra biologica foi sequenciada em multiplas lanes, repetir o mesmo nome na coluna `sample` -- o pipeline concatenara os reads automaticamente.
- Nao usar espacos nos nomes das amostras (usar underscores).
- Garantir que os arquivos FASTQ estejam compactados (`.fastq.gz`).

---

## Etapa 5 -- Obtencao do Genoma de Referencia

O pipeline precisa de um genoma de referencia (FASTA) e de uma anotacao (GTF) para realizar o alinhamento e a quantificacao.

### 5.1 Fontes comuns de genomas

| Organismo | Fonte | URL |
|-----------|-------|-----|
| Humano (GRCh38) | Ensembl | https://www.ensembl.org/Homo_sapiens/Info/Index |
| Camundongo (GRCm39) | Ensembl | https://www.ensembl.org/Mus_musculus/Info/Index |
| Arabidopsis | Ensembl Plants | https://plants.ensembl.org/Arabidopsis_thaliana/Info/Index |
| Outros | NCBI | https://www.ncbi.nlm.nih.gov/datasets/genome/ |

### 5.2 Exemplo: Baixar genoma humano (Ensembl)

```bash
cd data/reference/

# Baixar o genoma FASTA (primary assembly)
wget https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Baixar a anotacao GTF
wget https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz
gunzip Homo_sapiens.GRCh38.112.gtf.gz

cd ../..
```

### 5.3 Usar genomas pre-configurados do iGenomes (alternativa)

O pipeline suporta genomas pre-configurados do iGenomes atraves do parametro `--genome`. No entanto, essa abordagem **nao e recomendada** pela equipe do nf-core, pois os genomas podem estar desatualizados.

```bash
# Exemplo (nao recomendado)
--genome GRCh38
```

A abordagem recomendada e fornecer os arquivos explicitamente com `--fasta` e `--gtf`.

---

## Etapa 6 -- Execucao do Pipeline nf-core/rnaseq

### 6.1 Verificar se Docker esta em execucao

Antes de iniciar, confirmar que o Docker esta rodando:

```bash
docker info
```

Se houver erro de conexao, iniciar o Docker:
- **Linux:** `sudo systemctl start docker`
- **macOS/Windows:** abrir o aplicativo Docker Desktop.

### 6.2 (Opcional) Executar o teste do pipeline

Antes de rodar seus dados, e recomendavel executar um teste para garantir que tudo esta configurado corretamente:

```bash
nextflow run nf-core/rnaseq \
    -profile test,docker \
    --outdir results/test
```

Este teste usa um pequeno conjunto de dados incluido no pipeline. Se finalizar sem erros, o ambiente esta pronto.

### 6.3 Executar o pipeline com seus dados

```bash
nextflow run nf-core/rnaseq \
    --input samplesheet.csv \
    --outdir results \
    --fasta data/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --gtf data/reference/Homo_sapiens.GRCh38.112.gtf \
    -profile docker
```

### 6.4 Parametros importantes

| Parametro | Descricao | Exemplo |
|-----------|-----------|---------|
| `--input` | Caminho para o samplesheet CSV | `samplesheet.csv` |
| `--outdir` | Diretorio de saida dos resultados | `results` |
| `--fasta` | Genoma de referencia em formato FASTA | `data/reference/genome.fa` |
| `--gtf` | Anotacao do genoma em formato GTF | `data/reference/genes.gtf` |
| `-profile` | Perfil de execucao (usar `docker`) | `docker` |
| `--aligner` | Alinhador a ser utilizado | `star_salmon` (padrao), `star_rsem`, `hisat2` |
| `--pseudo_aligner` | Pseudo-alinhador (sem alinhamento ao genoma) | `salmon` |
| `--skip_trimming` | Pular a etapa de trimagem | `true` ou `false` |
| `--extra_trimgalore_args` | Argumentos extras para Trim Galore | `'--clip_r1 10'` |
| `--min_mapped_reads` | Percentual minimo de reads mapeados | `5` (padrao) |

### 6.5 Exemplos de execucao com diferentes alinhadores

#### Usando STAR + Salmon (padrao, recomendado)

```bash
nextflow run nf-core/rnaseq \
    --input samplesheet.csv \
    --outdir results \
    --fasta data/reference/genome.fa \
    --gtf data/reference/genes.gtf \
    --aligner star_salmon \
    -profile docker
```

#### Usando apenas Salmon (pseudo-alinhamento, mais rapido)

```bash
nextflow run nf-core/rnaseq \
    --input samplesheet.csv \
    --outdir results \
    --fasta data/reference/genome.fa \
    --gtf data/reference/genes.gtf \
    --pseudo_aligner salmon \
    --skip_alignment \
    -profile docker
```

#### Usando HISAT2

```bash
nextflow run nf-core/rnaseq \
    --input samplesheet.csv \
    --outdir results \
    --fasta data/reference/genome.fa \
    --gtf data/reference/genes.gtf \
    --aligner hisat2 \
    -profile docker
```

### 6.6 Salvar os parametros em um arquivo YAML (recomendado)

Para reprodutibilidade, salvar todos os parametros em um arquivo `params.yaml`:

```yaml
input: 'samplesheet.csv'
outdir: 'results'
fasta: 'data/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
gtf: 'data/reference/Homo_sapiens.GRCh38.112.gtf'
aligner: 'star_salmon'
```

E executar com:

```bash
nextflow run nf-core/rnaseq \
    -profile docker \
    -params-file params.yaml
```

### 6.7 Retomar uma execucao interrompida

Se o pipeline falhar no meio da execucao (queda de energia, falta de memoria, etc.), usar o parametro `-resume` para continuar de onde parou:

```bash
nextflow run nf-core/rnaseq \
    -profile docker \
    -params-file params.yaml \
    -resume
```

O Nextflow reutilizara os resultados ja calculados, economizando tempo.

### 6.8 Executar em segundo plano

Para execucoes longas, usar `screen` ou `nohup`:

```bash
# Opcao 1: com nohup
nohup nextflow run nf-core/rnaseq \
    -profile docker \
    -params-file params.yaml \
    > log_rnaseq.txt 2>&1 &

# Opcao 2: com screen
screen -S rnaseq
nextflow run nf-core/rnaseq \
    -profile docker \
    -params-file params.yaml
# Para desanexar: Ctrl+A, depois D
# Para reconectar: screen -r rnaseq
```

---

## Etapa 7 -- Interpretacao dos Resultados

Apos a execucao bem-sucedida, os resultados estarao no diretorio especificado em `--outdir`.

### 7.1 Estrutura dos resultados

```
results/
├── fastqc/                   # Relatorios de qualidade dos reads brutos
│   ├── amostra1_R1_fastqc.html
│   └── ...
├── trimgalore/               # Reads apos trimagem
│   ├── amostra1_R1_trimmed.fastq.gz
│   └── ...
├── star_salmon/              # Resultados do alinhamento e quantificacao
│   ├── amostra1/
│   │   ├── amostra1.markdup.sorted.bam   # Alinhamento (BAM)
│   │   └── amostra1.markdup.sorted.bam.bai
│   └── ...
├── star_salmon/
│   └── salmon.merged.gene_counts.tsv     # **Matriz de contagem de genes**
│   └── salmon.merged.gene_tpm.tsv        # Valores TPM
│   └── salmon.merged.gene_counts_length_scaled.tsv
├── multiqc/                  # Relatorio consolidado
│   └── multiqc_report.html   # Abrir no navegador
└── pipeline_info/            # Informacoes sobre a execucao
    └── execution_report.html
```

### 7.2 Arquivos mais importantes

| Arquivo | Descricao | Uso |
|---------|-----------|-----|
| `multiqc/multiqc_report.html` | Relatorio completo de QC | Avaliar qualidade geral das amostras |
| `star_salmon/salmon.merged.gene_counts.tsv` | Matriz de contagem de genes (raw counts) | Entrada para analise diferencial |
| `star_salmon/salmon.merged.gene_tpm.tsv` | Valores TPM normalizados | Visualizacao e comparacoes |
| `star_salmon/salmon.merged.gene_counts_length_scaled.tsv` | Contagens corrigidas por comprimento | Alternativa para analise diferencial |

### 7.3 Verificacao de qualidade

Abrir o relatorio MultiQC no navegador:

```bash
# Linux
xdg-open results/multiqc/multiqc_report.html

# macOS
open results/multiqc/multiqc_report.html
```

Pontos a verificar no relatorio:

- **FastQC:** qualidade dos reads (Phred score > 20 na maioria das posicoes).
- **Trimming:** percentual de reads retidos apos trimagem (idealmente > 90%).
- **Alinhamento:** taxa de mapeamento (idealmente > 70% para eucariotos).
- **Contagem de genes:** distribuicao das contagens entre amostras.
- **PCA / Correlacao:** amostras do mesmo grupo devem agrupar juntas.

---

## Etapa 8 -- Analise de Expressao Diferencial (DEGs)

Apos obter as matrizes de contagem do pipeline nf-core/rnaseq, existem duas abordagens para realizar a analise de expressao diferencial.

### Abordagem A: Usando o pipeline nf-core/differentialabundance

O nf-core oferece um pipeline dedicado para analise diferencial que aceita diretamente as saidas do nf-core/rnaseq.

#### A.1 Preparar o samplesheet para analise diferencial

Criar o arquivo `samplesheet_diff.csv` com informacoes sobre as condicoes experimentais:

```csv
sample,condition
CONTROLE_REP1,controle
CONTROLE_REP2,controle
CONTROLE_REP3,controle
TRATAMENTO_REP1,tratamento
TRATAMENTO_REP2,tratamento
TRATAMENTO_REP3,tratamento
```

#### A.2 Preparar o arquivo de contrastes

Criar o arquivo `contrasts.csv` especificando quais comparacoes realizar:

```csv
id,variable,reference,target
tratamento_vs_controle,condition,controle,tratamento
```

- `id`: nome do contraste (descritivo).
- `variable`: nome da coluna no samplesheet que contem os grupos.
- `reference`: grupo de referencia (denominador no fold change).
- `target`: grupo de interesse (numerador no fold change).

#### A.3 Executar o pipeline differentialabundance

```bash
nextflow run nf-core/differentialabundance \
    --input samplesheet_diff.csv \
    --contrasts contrasts.csv \
    --matrix results/star_salmon/salmon.merged.gene_counts_length_scaled.tsv \
    --gtf data/reference/Homo_sapiens.GRCh38.112.gtf \
    --outdir results/degs \
    -profile docker
```

### Abordagem B: Usando R com DESeq2 (manual)

Para maior controle sobre a analise, usar DESeq2 diretamente em R.

#### B.1 Instalar dependencias no R

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "tximport", "AnnotationDbi"))
install.packages(c("readr", "dplyr", "ggplot2", "pheatmap"))
```

#### B.2 Script completo de analise com DESeq2

Criar o arquivo `scripts/deseq2_analysis.R`:

```r
# ============================================================
# Analise de Expressao Diferencial com DESeq2
# Entrada: matrizes de contagem do nf-core/rnaseq
# ============================================================

library(DESeq2)
library(readr)
library(dplyr)
library(ggplot2)
library(pheatmap)

# --------------------------------------------------
# 1. Carregar a matriz de contagem
# --------------------------------------------------
counts_file <- "results/star_salmon/salmon.merged.gene_counts.tsv"
counts_raw <- read_tsv(counts_file)

# A primeira coluna e o gene_id e a segunda e gene_name
gene_info <- counts_raw[, c("gene_id", "gene_name")]

# Montar a matriz de contagem (somente colunas numericas)
count_matrix <- as.matrix(counts_raw[, -c(1, 2)])
rownames(count_matrix) <- counts_raw$gene_id

# Arredondar para inteiros (necessario para DESeq2)
count_matrix <- round(count_matrix)

# --------------------------------------------------
# 2. Definir as condicoes experimentais
# --------------------------------------------------
# Ajustar conforme os nomes das suas colunas no arquivo de contagem
sample_names <- colnames(count_matrix)

condition <- ifelse(grepl("CONTROLE", sample_names), "controle", "tratamento")

col_data <- data.frame(
    row.names = sample_names,
    condition = factor(condition, levels = c("controle", "tratamento"))
)

# --------------------------------------------------
# 3. Criar o objeto DESeq2 e executar a analise
# --------------------------------------------------
dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = col_data,
    design = ~ condition
)

# Filtrar genes com baixa contagem (manter genes com ao menos 10 reads no total)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Executar a analise diferencial
dds <- DESeq(dds)

# Extrair os resultados
res <- results(dds, contrast = c("condition", "tratamento", "controle"))
res <- res[order(res$padj), ]

# --------------------------------------------------
# 4. Filtrar DEGs significativos
# --------------------------------------------------
# Criterios padrao: |log2FoldChange| > 1 e padj < 0.05
degs <- as.data.frame(res) %>%
    filter(!is.na(padj)) %>%
    filter(padj < 0.05 & abs(log2FoldChange) > 1)

# Adicionar nomes dos genes
degs$gene_id <- rownames(degs)
degs <- merge(degs, gene_info, by = "gene_id", all.x = TRUE)

# Classificar como up ou down regulado
degs$regulacao <- ifelse(degs$log2FoldChange > 0, "UP", "DOWN")

cat("Total de DEGs encontrados:", nrow(degs), "\n")
cat("UP-regulados:", sum(degs$regulacao == "UP"), "\n")
cat("DOWN-regulados:", sum(degs$regulacao == "DOWN"), "\n")

# --------------------------------------------------
# 5. Salvar resultados
# --------------------------------------------------
dir.create("results/degs", showWarnings = FALSE)

# Tabela completa de resultados (todos os genes)
write.csv(as.data.frame(res), "results/degs/resultados_completos.csv")

# Tabela com apenas os DEGs significativos
write.csv(degs, "results/degs/degs_significativos.csv", row.names = FALSE)

# DEGs upregulados
write.csv(
    degs %>% filter(regulacao == "UP"),
    "results/degs/degs_upregulados.csv",
    row.names = FALSE
)

# DEGs downregulados
write.csv(
    degs %>% filter(regulacao == "DOWN"),
    "results/degs/degs_downregulados.csv",
    row.names = FALSE
)

# --------------------------------------------------
# 6. Visualizacoes
# --------------------------------------------------
# Volcano plot
pdf("results/degs/volcano_plot.pdf", width = 10, height = 8)
res_df <- as.data.frame(res) %>%
    mutate(
        significant = ifelse(
            !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1,
            ifelse(log2FoldChange > 0, "UP", "DOWN"),
            "NS"
        )
    )

ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
    geom_point(alpha = 0.5, size = 1.5) +
    scale_color_manual(values = c("UP" = "#D73027", "DOWN" = "#4575B4", "NS" = "grey60")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
    labs(
        title = "Volcano Plot",
        x = "log2 Fold Change",
        y = "-log10 p-value",
        color = "Regulacao"
    ) +
    theme_minimal()
dev.off()

# PCA
pdf("results/degs/pca_plot.pdf", width = 8, height = 6)
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "condition") +
    theme_minimal() +
    labs(title = "PCA - Amostras")
dev.off()

# Heatmap dos top 50 DEGs
if (nrow(degs) > 0) {
    pdf("results/degs/heatmap_top_degs.pdf", width = 10, height = 12)
    top_genes <- head(degs$gene_id, 50)
    mat <- assay(vsd)[top_genes, ]
    mat <- mat - rowMeans(mat)
    pheatmap(
        mat,
        annotation_col = as.data.frame(colData(dds)["condition"]),
        scale = "row",
        show_rownames = TRUE,
        main = "Top 50 DEGs"
    )
    dev.off()
}

cat("\nAnalise concluida. Resultados salvos em results/degs/\n")
```

#### B.3 Executar o script

```bash
Rscript scripts/deseq2_analysis.R
```

---

## Resolucao de Problemas

### Docker nao inicia

```bash
# Linux: verificar o servico
sudo systemctl status docker
sudo systemctl start docker

# macOS/Windows: abrir o aplicativo Docker Desktop manualmente
```

### Erro de permissao com Docker

```bash
# Verificar se o usuario esta no grupo docker
groups $USER

# Se nao estiver, adicionar
sudo usermod -aG docker $USER
# Fazer logout e login novamente
```

### Pipeline falha por falta de memoria

Aumentar a memoria disponivel para o Docker (Docker Desktop > Settings > Resources) ou usar o parametro `--max_memory`:

```bash
nextflow run nf-core/rnaseq \
    -profile docker \
    -params-file params.yaml \
    --max_memory '12.GB' \
    --max_cpus 4
```

### Erro "No space left on device"

```bash
# Limpar containers e imagens Docker nao utilizados
docker system prune -a

# Verificar o espaco em disco
df -h
```

### Pipeline interrompido, como retomar

```bash
# Usar -resume para continuar de onde parou
nextflow run nf-core/rnaseq \
    -profile docker \
    -params-file params.yaml \
    -resume
```

### Limpar arquivos temporarios apos a execucao

```bash
# O diretorio work/ pode ocupar muito espaco. Remover apos confirmar os resultados:
rm -rf work/

# Manter os logs do Nextflow para referencia:
# .nextflow.log
```

---

## Referencias

- **nf-core/rnaseq:** https://nf-co.re/rnaseq
- **nf-core/differentialabundance:** https://nf-co.re/differentialabundance
- **Nextflow:** https://www.nextflow.io/docs/latest/index.html
- **Docker:** https://docs.docker.com/get-started/
- **DESeq2:** Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology. 2014;15(12):550.
- **STAR:** Dobin A et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013;29(1):15-21.
- **Salmon:** Patro R et al. Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods. 2017;14:417-419.
- **Ewels PA et al.** The nf-core framework for community-curated bioinformatics pipelines. Nature Biotechnology. 2020;38:276-278.

---

## Licenca

Este protocolo e disponibilizado sob a licenca MIT.
