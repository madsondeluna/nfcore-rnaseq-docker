# Pipeline de RNA-Seq com nf-core/rnaseq e Docker

Guia completo para executar análise de RNA-Seq - desde a instalação do Docker até a extração de genes diferencialmente expressos (DEGs) - utilizando o pipeline [nf-core/rnaseq](https://nf-co.re/rnaseq) com containers Docker.

---

## Sumário

1. [Visão Geral](#visão-geral)
2. [RNA-Seq *de novo* vs. RNA-Seq com genoma de referência](#rna-seq-de-novo-vs-rna-seq-com-genoma-de-referência)
3. [Pré-requisitos](#pré-requisitos)
4. [Etapa 1 - Instalação e Configuração do Docker](#etapa-1---instalação-e-configuração-do-docker)
5. [Etapa 2 - Instalação do Nextflow](#etapa-2---instalação-do-nextflow)
6. [Etapa 3 - Organização dos Dados](#etapa-3---organização-dos-dados)
7. [Etapa 4 - Preparação do Samplesheet](#etapa-4---preparação-do-samplesheet)
8. [Etapa 5 - Obtenção do Genoma de Referência](#etapa-5---obtenção-do-genoma-de-referência)
9. [Etapa 6 - Execução do Pipeline nf-core/rnaseq](#etapa-6---execução-do-pipeline-nf-corernaseq)
10. [Etapa 7 - Interpretação dos Resultados](#etapa-7---interpretação-dos-resultados)
11. [Etapa 8 - Análise de Expressão Diferencial (DEGs)](#etapa-8---análise-de-expressão-diferencial-degs)
12. [Resolução de Problemas](#resolução-de-problemas)
13. [Referências](#referências)

---

## Visão Geral

Este protocolo utiliza o pipeline **nf-core/rnaseq** (versão 3.18.0 ou superior) para processar dados brutos de RNA-Seq (arquivos FASTQ) e gerar matrizes de contagem de genes. Posteriormente, utiliza o pipeline **nf-core/differentialabundance** para identificar genes diferencialmente expressos (DEGs).

### O que o pipeline nf-core/rnaseq faz

O pipeline nf-core/rnaseq é dividido em **5 estágios principais**, cada um contendo múltiplas ferramentas que são executadas automaticamente. A imagem abaixo (mapa no estilo metrô) ilustra todas as rotas possíveis dentro do pipeline, com diferentes cores representando os diferentes métodos de alinhamento/quantificação disponíveis.

> Documentação oficial completa: [nf-core/rnaseq](https://nf-co.re/rnaseq/3.22.2/docs/usage/)

![Mapa do pipeline nf-core/rnaseq no estilo metró, mostrando os 5 estágios de processamento](https://raw.githubusercontent.com/nf-core/rnaseq/3.22.2//docs/images/nf-core-rnaseq_metro_map_grey_animated.svg)

As cores no mapa representam os diferentes métodos disponíveis:

| Cor | Método | Comando |
|-----|--------|---------|
| Verde | STAR + Salmon (padrão) | `--aligner star_salmon` |
| Azul | STAR + RSEM | `--aligner star_rsem` |
| Amarelo | HISAT2 (sem quantificação integrada) | `--aligner hisat2` |
| Rosa | Salmon (pseudo-alinhamento) | `--pseudo_aligner salmon` |
| Roxo | Kallisto (pseudo-alinhamento) | `--pseudo_aligner kallisto` |

---

### Estágio 1 - Pré-processamento

O primeiro estágio prepara os reads brutos para o alinhamento. Todas as ferramentas deste estágio são executadas independentemente do alinhador escolhido.

| Etapa | Ferramenta | O que faz |
|-------|------------|-----------|
| Concatenação de reads | `cat` | Se uma amostra foi sequenciada em múltiplas lanes, os arquivos FASTQ são concatenados automaticamente |
| Controle de qualidade inicial | [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | Gera relatórios de qualidade dos reads brutos (distribuição de bases, qualidade Phred, conteúdo GC, sequências duplicadas) |
| Inferência de strandedness | [Salmon](https://salmon.readthedocs.io/) + [fq](https://github.com/stjude-rust-labs/fq) | Quando `strandedness` é definido como `auto`, o pipeline sub-amostra 1 milhão de reads e usa o Salmon para inferir automaticamente a orientação da biblioteca |
| Extração de UMI | [UMI-tools](https://umi-tools.readthedocs.io/) | (Opcional) Extrai identificadores moleculares únicos (UMIs) dos reads, se aplicável ao protocolo de sequenciamento |
| Trimagem de adaptadores | [Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) ou [fastp](https://github.com/OpenGene/fastp) | Remove adaptadores de sequenciamento e bases de baixa qualidade das extremidades dos reads |
| Remoção de rRNA | [SortMeRNA](https://github.com/sortmerna/sortmerna) | (Opcional) Remove reads de RNA ribossômico que podem contaminar a biblioteca |
| Separação de genomas | [BBSplit](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) | (Opcional) Separa reads provenientes de diferentes organismos em amostras contaminadas ou de xenoenxerto |
| Controle de qualidade pós-trimagem | [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | Gera novos relatórios de qualidade após a trimagem para confirmar a melhoria |

Para pular a trimagem (se os dados já foram trimados):

```bash
nextflow run nf-core/rnaseq \
    --input samplesheet.csv \
    --outdir results \
    --fasta genome.fa \
    --gtf genes.gtf \
    --skip_trimming \
    -profile docker
```

Para usar fastp em vez de Trim Galore:

```bash
nextflow run nf-core/rnaseq \
    --input samplesheet.csv \
    --outdir results \
    --fasta genome.fa \
    --gtf genes.gtf \
    --trimmer fastp \
    -profile docker
```

> Documentação dos parâmetros de pré-processamento: [nf-core/rnaseq - Read trimming](https://nf-co.re/rnaseq/3.22.2/docs/usage/#adapter-trimming-options)

---

### Estágio 2 - Alinhamento ao genoma e quantificação

Neste estágio, os reads limpos são alinhados ao genoma de referência. O método padrão é **STAR + Salmon**, mas existem alternativas. Cada rota aparece com uma cor diferente no mapa.

| Rota | Alinhador | Quantificador | Comando |
|------|-----------|---------------|---------|
| Padrão (verde) | [STAR](https://github.com/alexdobin/STAR) | [Salmon](https://salmon.readthedocs.io/) | `--aligner star_salmon` |
| Alternativa (azul) | [STAR](https://github.com/alexdobin/STAR) | [RSEM](https://github.com/deweylab/RSEM) | `--aligner star_rsem` |
| Alternativa (amarelo) | [HISAT2](http://daehwankimlab.github.io/hisat2/) | - | `--aligner hisat2` |

Após o alinhamento, o pipeline executa automaticamente:

- **UMI-tools dedup** - Remove duplicatas de PCR baseadas em UMIs (se aplicável)
- **SAMtools** (sort, index, stats) - Ordena, indexa e gera estatísticas dos arquivos BAM
- **picard MarkDuplicates** - Marca reads duplicados no alinhamento
- **Salmon quant** (modo de mapeamento) - Quantifica a expressão gênica a partir do alinhamento STAR

```bash
# Rota padrão: STAR + Salmon (recomendado para a maioria dos casos)
nextflow run nf-core/rnaseq \
    --input samplesheet.csv \
    --outdir results \
    --fasta genome.fa \
    --gtf genes.gtf \
    --aligner star_salmon \
    -profile docker
```

> Documentação das opções de alinhamento: [nf-core/rnaseq - Alignment options](https://nf-co.re/rnaseq/3.22.2/docs/usage/#alignment-options)

---

### Estágio 3 - Pseudo-alinhamento e quantificação

Alternativa mais rápida ao alinhamento tradicional. O pseudo-alinhamento **não gera arquivos BAM** - ele quantifica diretamente a expressão gênica a partir dos reads, sem mapeá-los posição a posição no genoma.

| Rota | Ferramenta | Comando |
|------|------------|---------|
| Rosa | [Salmon](https://salmon.readthedocs.io/) | `--pseudo_aligner salmon` |
| Roxo | [Kallisto](https://pachterlab.github.io/kallisto/) | `--pseudo_aligner kallisto` |

O pseudo-alinhamento é significativamente mais rápido e usa menos memória, sendo ideal para:
- Análises exploratórias iniciais
- Máquinas com recursos limitados
- Estudos focados apenas em quantificação (sem necessidade de BAM)

```bash
# Pseudo-alinhamento com Salmon (sem alinhamento ao genoma)
nextflow run nf-core/rnaseq \
    --input samplesheet.csv \
    --outdir results \
    --fasta genome.fa \
    --gtf genes.gtf \
    --pseudo_aligner salmon \
    --skip_alignment \
    -profile docker
```

> Documentação de quantificação: [nf-core/rnaseq - Quantification options](https://nf-co.re/rnaseq/3.22.2/docs/usage/#quantification-options)

---

### Estágio 4 - Pós-processamento

Após o alinhamento, o pipeline gera arquivos auxiliares e executa análises de qualidade do alinhamento.

| Etapa | Ferramenta | O que faz |
|-------|------------|-----------|
| Cobertura genômica | [BEDTools genomecov](https://bedtools.readthedocs.io/) | Calcula a cobertura de reads ao longo do genoma |
| Visualização de cobertura | [bedGraphToBigWig](https://genome.ucsc.edu/goldenPath/help/bigWig.html) | Converte arquivos de cobertura para o formato BigWig, visualizável em genome browsers (IGV, UCSC) |
| Montagem de transcritos | [StringTie](https://ccb.jhu.edu/software/stringtie/) | (Opcional) Monta transcritos a partir do alinhamento, útil para descoberta de novos transcritos |

Esses arquivos permitem visualizar a cobertura de sequenciamento em ferramentas como o [IGV (Integrative Genomics Viewer)](https://igv.org/).

---

### Estágio 5 - Controle de qualidade final

O último estágio avalia a qualidade geral do experimento de RNA-Seq a partir de múltiplas perspectivas.

| Etapa | Ferramenta | O que avalia |
|-------|------------|--------------|
| Relatório integrado | [MultiQC](https://multiqc.info/) | Consolida todas as métricas de QC em um único relatório HTML interativo |
| Complexidade da biblioteca | [Preseq](http://smithlabresearch.org/software/preseq/) | Estima a complexidade da biblioteca e prediz a cobertura com sequenciamento adicional |
| Qualidade do alinhamento | [RSeQC](http://rseqc.sourceforge.net/) (múltiplos módulos) | Avalia distribuição de reads por região gênica, saturação de junções, uniformidade de cobertura, inferência de strandedness |
| Qualidade por região | [Qualimap rnaseq](http://qualimap.conesalab.org/) | Avalia a qualidade do alinhamento por região genômica (exons, introns, intergênico) |
| Taxa de duplicação | [dupRadar](https://bioconductor.org/packages/dupRadar/) | Avalia se a taxa de duplicação é proporcional à expressão (esperado) ou indica viés de amplificação |
| PCA exploratório | [DESeq2](https://bioconductor.org/packages/DESeq2/) (apenas PCA) | Gera análise de componentes principais para visualizar agrupamento das amostras |
| Contaminação | [Kraken2](https://ccb.jhu.edu/software/kraken2/) / [Bracken](https://ccb.jhu.edu/software/bracken/) | (Opcional) Identifica possíveis contaminações por outros organismos nos reads não alinhados |

O relatório MultiQC é o ponto de partida para avaliar a qualidade do experimento. Para habilitá-lo com triagem de contaminação:

```bash
nextflow run nf-core/rnaseq \
    --input samplesheet.csv \
    --outdir results \
    --fasta genome.fa \
    --gtf genes.gtf \
    --aligner star_salmon \
    --contaminant_screening kraken2 \
    --kraken2_db /caminho/para/banco_kraken2 \
    -profile docker
```

> Documentação completa dos estágios: [nf-core/rnaseq - Output](https://nf-co.re/rnaseq/3.22.2/docs/output/)

---

### Fluxo geral simplificado

```mermaid
flowchart LR
    A["FASTQ brutos"] --> B["Controle de<br>Qualidade<br>FastQC"]
    B --> C["Trimagem<br>Trim Galore"]
    C --> D["Alinhamento<br>STAR / HISAT2"]
    D --> E["Quantificação<br>Salmon"]
    E --> F["Relatório<br>MultiQC"]
    F --> G["Matrizes de<br>Contagem"]
    G --> H["Análise Diferencial<br>DESeq2 / edgeR"]
    H --> I["Lista de DEGs"]

    style A fill:#4a6fa5,stroke:#2d4a7a,color:#ffffff
    style B fill:#6b8cae,stroke:#4a6fa5,color:#ffffff
    style C fill:#6b8cae,stroke:#4a6fa5,color:#ffffff
    style D fill:#6b8cae,stroke:#4a6fa5,color:#ffffff
    style E fill:#6b8cae,stroke:#4a6fa5,color:#ffffff
    style F fill:#6b8cae,stroke:#4a6fa5,color:#ffffff
    style G fill:#e8c547,stroke:#c4a432,color:#333333
    style H fill:#d4763a,stroke:#b35e28,color:#ffffff
    style I fill:#4caf50,stroke:#388e3c,color:#ffffff
```

---

## RNA-Seq *de novo* vs. RNA-Seq com genoma de referência

Antes de executar o pipeline, é fundamental entender as duas estratégias principais de análise de RNA-Seq, pois a escolha impacta diretamente nos resultados, na precisão e nos recursos computacionais necessários.

### Princípios de cada abordagem

#### RNA-Seq com genoma de referência (reference-based)

Nesta abordagem, os reads são **alinhados diretamente a um genoma de referência** já disponível para o organismo de estudo. As coordenadas do alinhamento, combinadas com a anotação do genoma (arquivo GTF/GFF), permitem quantificar a expressão de cada gene conhecido.

O fluxo segue a lógica:

```
Reads FASTQ --> Alinhamento ao genoma (STAR/HISAT2) --> Quantificação (Salmon/featureCounts) --> Contagens por gene
```

Esta é a abordagem utilizada pelo **nf-core/rnaseq** quando um genoma de referência é fornecido (como no caso de *Vigna unguiculata* a partir do Phytozome). É a estratégia adotada neste protocolo.

#### RNA-Seq *de novo* (de novo assembly)

Quando **não existe um genoma de referência** disponível para o organismo (ou quando o genoma existente é de baixa qualidade), os reads são montados *de novo* em transcritos completos sem a necessidade de um genoma. Ferramentas como **Trinity** realizam essa montagem a partir dos próprios reads.

O fluxo segue a lógica:

```
Reads FASTQ --> Montagem de novo (Trinity) --> Anotação dos transcritos (BLAST/Trinotate) --> Quantificação --> Contagens por transcrito
```

O nf-core não possui um pipeline dedicado para montagem *de novo* de transcriptomas. Para essa abordagem, é necessário utilizar ferramentas externas como [Trinity](https://github.com/trinityrnaseq/trinityrnaseq) e [Trinotate](https://github.com/Trinotate/Trinotate).

### Quando usar cada abordagem

| Cenário | Abordagem recomendada |
|---------|----------------------|
| Genoma de referência de boa qualidade disponível | Reference-based |
| Genoma de referência ausente ou de baixa qualidade | *De novo* |
| Interesse em descobrir novos transcritos/isoformas | *De novo* ou reference-based com StringTie |
| Análise de expressão diferencial (DEGs) | Reference-based (maior precisão) |
| Organismos não-modelo sem genoma publicado | *De novo* |
| Recursos computacionais limitados | Reference-based (menos custoso) |

### Comparação detalhada

| Critério | Reference-based | *De novo* |
|----------|----------------|-----------|
| **Requisito principal** | Genoma de referência + anotação (GTF/GFF) | Apenas reads FASTQ |
| **Ferramentas típicas de alinhamento/montagem** | STAR, HISAT2, Salmon | Trinity, rnaSPAdes, Trans-ABySS |
| **Precisão na quantificação** | Alta, pois utiliza coordenadas genômicas conhecidas | Moderada, depende da qualidade da montagem |
| **Detecção de novos transcritos** | Limitada (possível com StringTie) | Alta, pois monta transcritos sem viés de anotação prévia |
| **Tempo de processamento** | Moderado (horas) | Alto (dias, dependendo do volume de dados) |
| **Uso de memória RAM** | Moderado (16-32 GB típico) | Alto (64-128 GB ou mais para Trinity) |
| **Espaço em disco** | Moderado | Alto (montagens intermediárias são volumosas) |
| **Complexidade da análise** | Menor, pipeline bem definido | Maior, envolve montagem + anotação + filtragem |
| **Reprodutibilidade** | Alta, com genoma de referência fixo | Depende da versão do montador e dos parâmetros |
| **Sensibilidade a erros de sequenciamento** | Menor, o genoma serve como guia | Maior, erros podem gerar transcritos quiméricos |
| **Pipelines nf-core disponíveis** | nf-core/rnaseq | Não disponível (usar Trinity manualmente) |
| **Anotação funcional** | Já disponível no GTF/GFF | Necessita etapa adicional (Trinotate, BLAST, InterProScan) |

### Qual abordagem este protocolo utiliza

Este protocolo adota a abordagem **reference-based**, utilizando o genoma de *Vigna unguiculata* v1.2 (IT97K-499-35) disponível no [Phytozome v14](https://phytozome-next.jgi.doe.gov/info/Vunguiculata_v1_2) como referência. Esta é a escolha adequada porque:

- O genoma de *V. unguiculata* está bem montado e anotado
- Permite quantificação precisa por gene utilizando coordenadas genômicas
- É compatível com o pipeline nf-core/rnaseq
- Exige menos recursos computacionais do que a montagem *de novo*
- Facilita a análise de expressão diferencial downstream com DESeq2 ou edgeR

---

## Pré-requisitos

| Requisito | Mínimo | Recomendado |
|-----------|--------|-------------|
| RAM | 8 GB | 16 GB ou mais |
| Disco | 50 GB livres | 200 GB ou mais |
| CPU | 4 núcleos | 8 núcleos ou mais |
| Sistema Operacional | Linux, macOS ou Windows (via WSL2) | Linux |
| Java | versão 11 ou superior | versão 17 |
| Conexão com internet | necessária para download de containers e genomas | - |

---

## Etapa 1 - Instalação e Configuração do Docker

O Docker permite executar softwares em containers isolados, garantindo que todas as dependências estejam corretas e que o ambiente seja reprodutivel.

### 1.1 Instalação no Linux (Ubuntu/Debian)

```bash
# Atualizar os pacotes do sistema
sudo apt-get update

# Instalar dependências necessarias
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

### 1.2 Instalação no macOS

1. Acessar o site oficial: https://www.docker.com/products/docker-desktop/
2. Baixar o **Docker Desktop** para macOS (Apple Silicon ou Intel, conforme seu hardware).
3. Abrir o arquivo `.dmg` baixado.
4. Arrastar o ícone do Docker para a pasta **Applications**.
5. Abrir o Docker Desktop a partir do Launchpad ou da pasta Applications.
6. Aguardar a inicialização completa (o ícone da baleia na barra de menu ficará estável).

### 1.3 Instalação no Windows

1. Acessar https://www.docker.com/products/docker-desktop/
2. Baixar o **Docker Desktop** para Windows.
3. Executar o instalador e seguir as instruções na tela.
4. **Importante:** habilitar o WSL2 quando solicitado durante a instalação.
5. Reiniciar o computador se solicitado.
6. Abrir o Docker Desktop e aguardar a inicialização.

### 1.4 Configuração pós-instalação (Linux)

Por padrão no Linux, o Docker exige permissões de superusuário. Para executar sem `sudo`:

```bash
# Adicionar seu usuário ao grupo docker
sudo usermod -aG docker $USER

# Aplicar a mudança (necessario logout/login, ou executar)
newgrp docker
```

### 1.5 Verificar se o Docker está funcionando

```bash
# Verificar a versão instalada
docker --version

# Executar o container de teste
docker run hello-world
```

Se a mensagem "Hello from Docker!" aparecer, a instalação foi bem-sucedida.

### 1.6 Configurar recursos do Docker

Para análises de RNA-Seq, é recomendável alocar recursos adequados ao Docker:

- **Docker Desktop (macOS/Windows):** abrir Docker Desktop > Settings > Resources > ajustar CPU, memória e disco conforme necessidade.
- **Linux:** o Docker utiliza os recursos do host diretamente, sem necessidade de configuração adicional.

---

## Etapa 2 - Instalação do Nextflow

O Nextflow é o motor de workflow que executa o pipeline nf-core/rnaseq. Ele orquestra os containers Docker automaticamente.

### 2.1 Verificar se o Javaestá instalado

```bash
java -version
```

Se o Java não estiver instalado ou for inferior à versão 11:

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

# Tornar o binario executável
chmod +x nextflow

# Mover para um diretório no PATH do sistema
sudo mv nextflow /usr/local/bin/
```

### 2.3 Verificar a instalação

```bash
nextflow -version
```

### 2.4 (Opcional) Instalar as ferramentas de linha de comando do nf-core

As ferramentas nf-core não são obrigatórias para rodar o pipeline, mas facilitam tarefas como baixar pipelines e gerar samplesheets.

```bash
pip install nf-core
```

---

## Etapa 3 - Organização dos Dados

### 3.1 Estrutura de diretórios recomendada

Criar a seguinte estrutura antes de iniciar a análise:

```
projeto_rnaseq/
├── data/
│   ├── raw_fastq/          # Arquivos FASTQ brutos aqui
│   │   ├── amostra1_R1.fastq.gz
│   │   ├── amostra1_R2.fastq.gz
│   │   ├── amostra2_R1.fastq.gz
│   │   ├── amostra2_R2.fastq.gz
│   │   └── ...
│   └── reference/           # Genoma de referência
│       ├── genome.fa
│       └── genes.gtf
├── samplesheet.csv           # Samplesheet de entrada
├── results/                  # Resultados do pipeline
└── work/                     # Diretório de trabalho do Nextflow
```

### 3.2 Criar os diretórios

```bash
mkdir -p projeto_rnaseq/{data/{raw_fastq,reference},results}
cd projeto_rnaseq
```

### 3.3 Transferir os arquivos FASTQ

Copiar os arquivos FASTQ brutos para o diretório `data/raw_fastq/`. Os arquivos devem estar compactados com gzip (extensao `.fastq.gz` ou `.fq.gz`).

```bash
# Exemplo: copiar de um HD externo
cp /media/hd_externo/fastq/*.fastq.gz data/raw_fastq/

# Exemplo: baixar de um servidor via scp
scp usuário@servidor:/caminho/dos/dados/*.fastq.gz data/raw_fastq/
```

---

## Etapa 4 - Preparação do Samplesheet

O samplesheet é um arquivo CSV que informa ao pipeline quais amostras processar, onde estão os arquivos FASTQ e qual a orientação da biblioteca (strandedness).

### 4.1 Formato do samplesheet

O arquivo deve conter exatamente 4 colunas, separadas por vírgula, com cabeçalho:

| Coluna | Descrição |
|--------|-----------|
| `sample` | Nome da amostra (identificador único por amostra biologica) |
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
| `unstranded` | Quando a biblioteca não preserva informação de fita (ex.: TruSeq unstranded) |
| `forward` | Biblioteca strand-specific no sentido forward (ex.: Ligation, dUTP - SMARTer) |
| `reverse` | Biblioteca strand-specific no sentido reverse (ex.: TruSeq Stranded, maioria dos kits Illumina atuais) |

### 4.5 Regras importantes

- Usar **caminhos absolutos** (completos) para os arquivos FASTQ.
- Se uma mesma amostra biológica foi sequenciada em múltiplas lanes, repetir o mesmo nome na coluna `sample` - o pipeline concatenará os reads automaticamente.
- Não usar espacos nos nomes das amostras (usar underscores).
- Garantir que os arquivos FASTQ estejam compactados (`.fastq.gz`).

---

## Etapa 5 - Obtenção do Genoma de Referência

O pipeline precisa de um genoma de referência (FASTA) e de uma anotação (GTF/GFF) para realizar o alinhamento e a quantificação. Neste protocolo, utilizamos o genoma de **Vigna unguiculata** (feijão-caupi), disponível no [Phytozome v14](https://phytozome-next.jgi.doe.gov/).

### 5.1 Fonte do genoma

| Organismo | Versão | Fonte | URL |
|-----------|--------|-------|-----|
| *Vigna unguiculata* (feijão-caupi) | v1.2 (IT97K-499-35) | Phytozome v14 | https://phytozome-next.jgi.doe.gov/info/Vunguiculata_v1_2 |

### 5.2 Como baixar o genoma do Phytozome

O Phytozome exige login (gratuito) para download dos dados. Siga os passos abaixo:

1. Criar uma conta gratuita em https://phytozome-next.jgi.doe.gov/ (caso ainda não possua).
2. Acessar a página de *Vigna unguiculata*: https://phytozome-next.jgi.doe.gov/info/Vunguiculata_v1_2
3. Clicar em **Download** e baixar os seguintes arquivos:

| Arquivo | Descrição | Nome esperado |
|---------|-----------|---------------|
| Genoma FASTA | Sequências dos cromossomos/scaffolds | `Vunguiculata_469_v1.0.fa.gz` |
| Anotação GFF3 | Coordenadas dos genes e transcritos | `Vunguiculata_469_v1.2.gene.gff3.gz` |

4. Mover os arquivos para o diretório de referência e descompactar:

```bash
cd data/reference/

# Descompactar o genoma FASTA
gunzip Vunguiculata_469_v1.0.fa.gz

# Descompactar a anotação GFF3
gunzip Vunguiculata_469_v1.2.gene.gff3.gz

cd ../..
```

> **Nota:** O Phytozome fornece anotação no formato GFF3. O pipeline nf-core/rnaseq aceita tanto GTF quanto GFF3 através do parâmetro `--gff` (em vez de `--gtf`). Se preferir converter para GTF, utilize o [gffread](https://github.com/gpertea/gffread):
>
> ```bash
> gffread Vunguiculata_469_v1.2.gene.gff3 -T -o Vunguiculata_469_v1.2.gtf
> ```

### 5.3 Usar genomas pré-configurados do iGenomes (alternativa)

O pipeline suporta genomas pré-configurados do iGenomes através do parâmetro `--genome`. No entanto, *Vigna unguiculata* **não está disponível** no iGenomes, portanto é necessário fornecer os arquivos manualmente com `--fasta` e `--gff` (ou `--gtf`).

---

## Etapa 6 - Execução do Pipeline nf-core/rnaseq

### 6.1 Verificar se Docker está em execução

Antes de iniciar, confirmar que o Docker está rodando:

```bash
docker info
```

Se houver erro de conexão, iniciar o Docker:
- **Linux:** `sudo systemctl start docker`
- **macOS/Windows:** abrir o aplicativo Docker Desktop.

### 6.2 (Opcional) Executar o teste do pipeline

Antes de rodar seus dados, é recomendável executar um teste para garantir que tudo está configurado corretamente:

```bash
nextflow run nf-core/rnaseq \
    -profile test,docker \
    --outdir results/test
```

Este teste usa um pequeno conjunto de dados incluído no pipeline. Se finalizar sem erros, o ambiente está pronto.

### 6.3 Executar o pipeline com seus dados

```bash
nextflow run nf-core/rnaseq \
    --input samplesheet.csv \
    --outdir results \
    --fasta data/reference/Vunguiculata_469_v1.0.fa \
    --gff data/reference/Vunguiculata_469_v1.2.gene.gff3 \
    -profile docker
```

### 6.4 Parâmetros importantes

| Parâmetro | Descrição | Exemplo |
|-----------|-----------|---------|
| `--input` | Caminho para o samplesheet CSV | `samplesheet.csv` |
| `--outdir` | Diretório de saída dos resultados | `results` |
| `--fasta` | Genoma de referência em formato FASTA | `data/reference/genome.fa` |
| `--gtf` | Anotacao do genoma em formato GTF | `data/reference/genes.gtf` |
| `-profile` | Perfil de execução (usar `docker`) | `docker` |
| `--aligner` | Alinhador a ser utilizado | `star_salmon` (padrão), `star_rsem`, `hisat2` |
| `--pseudo_aligner` | Pseudo-alinhador (sem alinhamento ao genoma) | `salmon` |
| `--skip_trimming` | Pular a etapa de trimagem | `true` ou `false` |
| `--extra_trimgalore_args` | Argumentos extras para Trim Galore | `'--clip_r1 10'` |
| `--min_mapped_reads` | Percentual mínimo de reads mapeados | `5` (padrão) |

### 6.5 Exemplos de execução com diferentes alinhadores

#### Usando STAR + Salmon (padrão, recomendado)

```bash
nextflow run nf-core/rnaseq \
    --input samplesheet.csv \
    --outdir results \
    --fasta data/reference/genome.fa \
    --gtf data/reference/genes.gtf \
    --aligner star_salmon \
    -profile docker
```

#### Usando apenas Salmon (pseudo-alinhamento, mais rápido)

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

### 6.6 Salvar os parâmetros em um arquivo YAML (recomendado)

Para reprodutibilidade, salvar todos os parâmetros em um arquivo `params.yaml`:

```yaml
input: 'samplesheet.csv'
outdir: 'results'
fasta: 'data/reference/Vunguiculata_469_v1.0.fa'
gff: 'data/reference/Vunguiculata_469_v1.2.gene.gff3'
aligner: 'star_salmon'
```

E executar com:

```bash
nextflow run nf-core/rnaseq \
    -profile docker \
    -params-file params.yaml
```

### 6.7 Retomar uma execução interrompida

Se o pipeline falhar no meio da execução (queda de energia, falta de memória, etc.), usar o parâmetro `-resume` para continuar de onde parou:

```bash
nextflow run nf-core/rnaseq \
    -profile docker \
    -params-file params.yaml \
    -resume
```

O Nextflow reutilizara os resultados já calculados, economizando tempo.

### 6.8 Executar em segundo plano

Para execucoes longas, usar `screen` ou `nohup`:

```bash
# Opção 1: com nohup
nohup nextflow run nf-core/rnaseq \
    -profile docker \
    -params-file params.yaml \
    > log_rnaseq.txt 2>&1 &

# Opção 2: com screen
screen -S rnaseq
nextflow run nf-core/rnaseq \
    -profile docker \
    -params-file params.yaml
# Para desanexar: Ctrl+A, depois D
# Para reconectar: screen -r rnaseq
```

---

## Etapa 7 - Interpretação dos Resultados

Após a execução bem-sucedida, os resultados estarão no diretório especificado em `--outdir`.

### 7.1 Estrutura dos resultados

```
results/
├── fastqc/                   # Relatorios de qualidade dos reads brutos
│   ├── amostra1_R1_fastqc.html
│   └── ...
├── trimgalore/               # Reads após trimagem
│   ├── amostra1_R1_trimmed.fastq.gz
│   └── ...
├── star_salmon/              # Resultados do alinhamento e quantificação
│   ├── amostra1/
│   │   ├── amostra1.markdup.sorted.bam   # Alinhamento (BAM)
│   │   └── amostra1.markdup.sorted.bam.bai
│   └── ...
├── star_salmon/
│   └── salmon.merged.gene_counts.tsv     # **Matriz de contagem de genes**
│   └── salmon.merged.gene_tpm.tsv        # Valores TPM
│   └── salmon.merged.gene_counts_length_scaled.tsv
├── multiqc/                  # Relatório consolidado
│   └── multiqc_report.html   # Abrir no navegador
└── pipeline_info/            # Informacoes sobre a execução
    └── execution_report.html
```

### 7.2 Arquivos mais importantes

| Arquivo | Descrição | Uso |
|---------|-----------|-----|
| `multiqc/multiqc_report.html` | Relatório completo de QC | Avaliar qualidade geral das amostras |
| `star_salmon/salmon.merged.gene_counts.tsv` | Matriz de contagem de genes (raw counts) | Entrada para análise diferencial |
| `star_salmon/salmon.merged.gene_tpm.tsv` | Valores TPM normalizados | Visualizacao e comparações |
| `star_salmon/salmon.merged.gene_counts_length_scaled.tsv` | Contagens corrigidas por comprimento | Alternativa para análise diferencial |

### 7.3 Verificacao de qualidade

Abrir o relatório MultiQC no navegador:

```bash
# Linux
xdg-open results/multiqc/multiqc_report.html

# macOS
open results/multiqc/multiqc_report.html
```

Pontos a verificar no relatorio:

- **FastQC:** qualidade dos reads (Phred score > 20 na maioria das posições).
- **Trimming:** percentual de reads retidos após trimagem (idealmente > 90%).
- **Alinhamento:** taxa de mapeamento (idealmente > 70% para eucariotos).
- **Contagem de genes:** distribuicao das contagens entre amostras.
- **PCA / Correlacao:** amostras do mesmo grupo devem agrupar juntas.

---

## Etapa 8 - Análise de Expressão Diferencial (DEGs)

Após obter as matrizes de contagem do pipeline nf-core/rnaseq, existem duas abordagens para realizar a análise de expressão diferencial.

### Abordagem A: Usando o pipeline nf-core/differentialabundance

O nf-core oferece um pipeline dedicado para análise diferencial que aceita diretamente as saidas do nf-core/rnaseq.

#### A.1 Preparar o samplesheet para análise diferencial

Criar o arquivo `samplesheet_diff.csv` com informações sobre as condições experimentais:

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

Criar o arquivo `contrasts.csv` especificando quais comparações realizar:

```csv
id,variable,reference,target
tratamento_vs_controle,condition,controle,tratamento
```

- `id`: nome do contraste (descritivo).
- `variable`: nome da coluna no samplesheet que contém os grupos.
- `reference`: grupo de referência (denominador no fold change).
- `target`: grupo de interesse (numerador no fold change).

#### A.3 Executar o pipeline differentialabundance

```bash
nextflow run nf-core/differentialabundance \
    --input samplesheet_diff.csv \
    --contrasts contrasts.csv \
    --matrix results/star_salmon/salmon.merged.gene_counts_length_scaled.tsv \
    --gff data/reference/Vunguiculata_469_v1.2.gene.gff3 \
    --outdir results/degs \
    -profile docker
```

### Abordagem B: Usando R com DESeq2 (manual)

Para maior controle sobre a análise, usar DESeq2 diretamente em R.

#### B.1 Instalar dependências no R

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "tximport", "AnnotationDbi"))
install.packages(c("readr", "dplyr", "ggplot2", "pheatmap"))
```

#### B.2 Script completo de análise com DESeq2

Criar o arquivo `scripts/deseq2_analysis.R`:

```r
# ============================================================
# Análise de Expressão Diferencial com DESeq2
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

# A primeira coluna é o gene_id e a segunda é gene_name
gene_info <- counts_raw[, c("gene_id", "gene_name")]

# Montar a matriz de contagem (somente colunas numéricas)
count_matrix <- as.matrix(counts_raw[, -c(1, 2)])
rownames(count_matrix) <- counts_raw$gene_id

# Arredondar para inteiros (necessario para DESeq2)
count_matrix <- round(count_matrix)

# --------------------------------------------------
# 2. Definir as condições experimentais
# --------------------------------------------------
# Ajustar conforme os nomes das suas colunas no arquivo de contagem
sample_names <- colnames(count_matrix)

condition <- ifelse(grepl("CONTROLE", sample_names), "controle", "tratamento")

col_data <- data.frame(
    row.names = sample_names,
    condition = factor(condition, levels = c("controle", "tratamento"))
)

# --------------------------------------------------
# 3. Criar o objeto DESeq2 e executar a análise
# --------------------------------------------------
dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = col_data,
    design = ~ condition
)

# Filtrar genes com baixa contagem (manter genes com ao menos 10 reads no total)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Executar a análise diferencial
dds <- DESeq(dds)

# Extrair os resultados
res <- results(dds, contrast = c("condition", "tratamento", "controle"))
res <- res[order(res$padj), ]

# --------------------------------------------------
# 4. Filtrar DEGs significativos
# --------------------------------------------------
# Critérios padrão: |log2FoldChange| > 1 e padj < 0.05
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
# 6. Visualizações
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

cat("\nAnálise concluída. Resultados salvos em results/degs/\n")
```

#### B.3 Executar o script

```bash
Rscript scripts/deseq2_analysis.R
```

---

## Resolução de Problemas

### Docker não inicia

```bash
# Linux: verificar o serviço
sudo systemctl status docker
sudo systemctl start docker

# macOS/Windows: abrir o aplicativo Docker Desktop manualmente
```

### Erro de permissão com Docker

```bash
# Verificar se o usuárioestá no grupo docker
groups $USER

# Se não estiver, adicionar
sudo usermod -aG docker $USER
# Fazer logout e login novamente
```

### Pipeline falha por falta de memória

Aumentar a memória disponível para o Docker (Docker Desktop > Settings > Resources) ou usar o parâmetro `--max_memory`:

```bash
nextflow run nf-core/rnaseq \
    -profile docker \
    -params-file params.yaml \
    --max_memory '12.GB' \
    --max_cpus 4
```

### Erro "No space left on device"

```bash
# Limpar containers e imagens Docker não utilizados
docker system prune -a

# Verificar o espaço em disco
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

### Limpar arquivos temporários após a execução

```bash
# O diretório work/ pode ocupar muito espaco. Remover após confirmar os resultados:
rm -rf work/

# Manter os logs do Nextflow para referência:
# .nextflow.log
```

---

## Referências

- **nf-core/rnaseq:** https://nf-co.re/rnaseq
- **nf-core/differentialabundance:** https://nf-co.re/differentialabundance
- **Nextflow:** https://www.nextflow.io/docs/latest/index.html
- **Docker:** https://docs.docker.com/get-started/
- **DESeq2:** Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology. 2014;15(12):550.
- **STAR:** Dobin A et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013;29(1):15-21.
- **Salmon:** Patro R et al. Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods. 2017;14:417-419.
- **Ewels PA et al.** The nf-core framework for community-curated bioinformatics pipelines. Nature Biotechnology. 2020;38:276-278.

---

## Licença

Este protocolo é disponibilizado sob a licença MIT.
