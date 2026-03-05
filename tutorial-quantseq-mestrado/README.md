# Tutorial Prático: RNA-Seq com nf-core/rnaseq e Docker

Tutorial completo para alunos de pós-graduação que estão aprendendo a montar e anotar dados de RNA-Seq do zero, utilizando o pipeline [nf-core/rnaseq](https://nf-co.re/rnaseq) dentro de containers Docker.

> Este tutorial é baseado no artigo: [QuantSeq RNAseq com nf-core (Thomas Singer, 2023)](https://tomsing1.github.io/blog/posts/nextflow-core-quantseq-1-settings/)

**Sistema utilizado neste tutorial:** Ubuntu Linux, processador Intel (x86_64)

---

## O que você vai aprender

1. O que é Docker e por que ele é essencial para bioinformática reprodutível
2. Como instalar o Docker Engine no Ubuntu
3. Como instalar o Nextflow
4. Como verificar se tudo está funcionando
5. Como baixar os dados de sequenciamento públicos (FASTQs do SRA)
6. Como montar o samplesheet e configurar o pipeline
7. Como executar o pipeline nf-core/rnaseq com Docker
8. Como interpretar os resultados

---

## Sumário

1. [Contexto: o que é QuantSeq RNA-Seq](#1-contexto-o-que-é-quantseq-rna-seq)
2. [Por que usar Docker em bioinformática](#2-por-que-usar-docker-em-bioinformática)
3. [O que você precisa instalar -e o que não precisa](#3-o-que-você-precisa-instalar--e-o-que-não-precisa)
4. [Parte 1 - Instalando o Docker Engine no Ubuntu](#parte-1---instalando-o-docker-engine-no-ubuntu)
5. [Parte 2 - Instalando o Nextflow](#parte-2---instalando-o-nextflow)
6. [Parte 3 - Verificando se tudo está funcionando](#parte-3---verificando-se-tudo-está-funcionando)
7. [Parte 4 - Baixando os dados do tutorial (SRA)](#parte-4---baixando-os-dados-do-tutorial-sra)
8. [Parte 5 - Preparando o samplesheet](#parte-5---preparando-o-samplesheet)
9. [Parte 6 - Executando o pipeline nf-core/rnaseq](#parte-6---executando-o-pipeline-nf-corernaseq)
10. [Parte 7 - Entendendo os resultados](#parte-7---entendendo-os-resultados)
11. [Solução de problemas comuns](#solução-de-problemas-comuns)

---

## 1. Contexto: o que é QuantSeq RNA-Seq

O **QuantSeq 3' mRNA-Seq** (da empresa Lexogen) é um protocolo de sequenciamento de RNA que captura apenas a região 3' dos transcritos. Isso reduz o custo e o tempo de processamento em comparação ao RNA-Seq convencional (que captura o transcrito inteiro), mantendo a capacidade de medir expressão gênica. Vale lembrar que o QuantSeq é diferente do RNA-Seq tradicional, e isso tem implicações para o processamento dos dados. Estamos usando esses dados devido ao seu tamanho reduzido e características específicas, o que torna o tutorial mais acessível.

**Características importantes:**
- Dados **single-end** (apenas uma direção de leitura, sem R2)
- Biblioteca **strand-specific** (sentido forward)
- Reads mais curtos (50-100 bp) concentrados na extremidade 3' dos genes
- Pode conter **UMIs** (Unique Molecular Identifiers) em alguns protocolos

**Dados que vamos usar neste tutorial:**

O tutorial original usa dois conjuntos de dados públicos de camundongo disponíveis no NCBI SRA:

| Dataset | Artigo | Acesso SRA | Acesso GEO | Detalhes |
|---------|--------|-----------|-----------|----------|
| Xia et al, 2021 | Sem UMIs | SRP282921 | GSE158152 | 36 amostras, camundongo |
| Nugent et al, 2020 | Com UMIs | SRP213880 | GSE134031 | Camundongo |

Nós vamos baixar e processar o dataset **Xia et al (SRP282921)**, que é o mais simples (sem UMIs).

---

## 2. Por que usar Docker em bioinformática

### O problema que o Docker resolve

Imagine que você instala um programa de bioinformática hoje e ele funciona perfeitamente. Daqui a um ano, um colega tenta reproduzir sua análise no computador dele -mas o programa não funciona, porque a versão do Python é diferente, ou uma biblioteca foi atualizada, ou o sistema operacional é outro.

Esse problema de **reprodutibilidade** é um dos maiores desafios da bioinformática moderna.

### O que é um container Docker

O Docker resolve isso criando **containers**: pacotes completos e isolados que contêm o programa, todas as suas dependências, as versões exatas das bibliotecas e até o sistema operacional mínimo necessário para rodar. É como uma "caixa fechada" que funciona identicamente em qualquer computador que tenha Docker instalado.

![Diagrama Docker](tutorial-quantseq-mestrado/docker-conceito.svg)

### Por que o nf-core usa Docker

O pipeline nf-core/rnaseq usa **dezenas de ferramentas diferentes** (STAR, Salmon, FastQC, Trim Galore, MultiQC, etc.). Em vez de você instalar cada uma manualmente (o que levaria horas e geraria conflitos de versão), o Nextflow baixa automaticamente os containers Docker de cada ferramenta quando necessário. 

Você só precisa ter o Docker instalado. O resto é automático.

---

## 3. O que você precisa instalar -e o que não precisa

Esta é uma das maiores vantagens de usar Docker com nf-core: a lista do que você realmente precisa instalar na sua máquina é mínima.

### Precisa instalar (3 coisas apenas):

| Software | Por quê precisa estar na máquina |
|----------|----------------------------------|
| **Docker Engine** | É o motor que cria e executa os containers |
| **Java** | O Nextflow é escrito em Java e precisa dele para iniciar |
| **Nextflow** | É quem lê o pipeline e orquestra a execução dos containers |

### NÃO precisa instalar:

- STAR, Salmon, HISAT2 -rodam dentro de containers Docker
- FastQC, Trim Galore, MultiQC -idem
- Python, R, Perl -as versões corretas já estão dentro dos containers
- nf-core tools -opcional, usaremos apenas para baixar os dados

O Nextflow baixa cada container automaticamente na primeira vez que o processo correspondente é executado e os armazena em cache para execuções futuras.

---

## Parte 1 - Instalando o Docker Engine no Ubuntu

O **Docker Engine** é a versão de linha de comando do Docker, ideal para Linux. Não é necessário instalar o Docker Desktop (que é uma interface gráfica voltada para macOS e Windows).

### 1.1 Remover versões antigas (se houver)

```bash
sudo apt-get remove docker docker-engine docker.io containerd runc
```

Esse comando remove instalações anteriores do Docker para evitar conflitos. Não se preocupe se aparecer a mensagem "package not found" -significa que não havia versão antiga instalada.

### 1.2 Instalar o Docker Engine

```bash
# Atualizar a lista de pacotes
sudo apt-get update

# Instalar dependências necessárias para o próximo passo
sudo apt-get install -y \
    ca-certificates \
    curl \
    gnupg \
    lsb-release
```

Esses pacotes permitem que o sistema baixe software de fontes HTTPS e verifique assinaturas digitais.

```bash
# Adicionar a chave GPG oficial do Docker
# Isso garante que o software baixado vem realmente do Docker Inc.
sudo mkdir -p /etc/apt/keyrings
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | \
    sudo gpg --dearmor -o /etc/apt/keyrings/docker.gpg
```

```bash
# Registrar o repositório oficial do Docker como fonte de pacotes
echo \
  "deb [arch=amd64 signed-by=/etc/apt/keyrings/docker.gpg] \
  https://download.docker.com/linux/ubuntu \
  $(lsb_release -cs) stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
```

Note que usamos `arch=amd64` porque o processador é Intel (x86_64). Se fosse ARM, seria `arch=arm64`.

```bash
# Instalar o Docker Engine
sudo apt-get update
sudo apt-get install -y docker-ce docker-ce-cli containerd.io docker-compose-plugin
```

### 1.3 Permitir usar Docker sem sudo

Por padrão no Linux, o Docker exige permissões de administrador (`sudo`) a cada comando. Para simplificar:

```bash
# Adicionar seu usuário ao grupo docker
sudo usermod -aG docker $USER

# Aplicar a mudança na sessão atual sem precisar fazer logout
newgrp docker
```

**Por que isso é necessário:** o Docker comunica com um processo em segundo plano (o "daemon") através de um socket que por padrão só aceita conexões do root. Ao adicionar seu usuário ao grupo `docker`, você ganha permissão para se comunicar com esse socket diretamente.

### 1.4 Iniciar o Docker e configurar para iniciar automaticamente

```bash
# Iniciar o Docker agora
sudo systemctl start docker

# Configurar para iniciar automaticamente ao ligar o computador
sudo systemctl enable docker

# Verificar se está rodando
sudo systemctl status docker
```

Você deve ver `Active: active (running)` na saída do último comando.

---

## Parte 2 - Instalando o Nextflow

O **Nextflow** é o motor de workflow que orquestra a execução do pipeline. Ele lê as instruções do nf-core/rnaseq, baixa os containers Docker de cada ferramenta e executa tudo na ordem correta.

O Nextflow em si roda diretamente na sua máquina (não dentro de um container), por isso precisa ser instalado localmente junto com o Java.

### 2.1 Instalar o Java

```bash
# Instalar o Java Development Kit (versão 17 recomendada)
sudo apt-get install -y default-jdk

# Verificar a instalação
java -version
```

Você deve ver algo como:
```
openjdk version "17.0.x" ...
```

### 2.2 Instalar o Nextflow

```bash
# Baixar o instalador do Nextflow
curl -s https://get.nextflow.io | bash

# O comando acima cria um arquivo executável "nextflow" no diretório atual
# Mover para /usr/local/bin/ para poder usá-lo de qualquer pasta
chmod +x nextflow
sudo mv nextflow /usr/local/bin/

# Verificar a instalação
nextflow -version
```

Você deve ver algo como:
```
      N E X T F L O W
      version 24.x.x build ...
      http://nextflow.io
```

**Por que mover para `/usr/local/bin/`:** esse diretório está incluído no PATH do sistema, o que permite chamar `nextflow` de qualquer pasta sem precisar escrever o caminho completo.

---

## Parte 3 - Verificando se tudo está funcionando

Antes de rodar seus dados reais, execute estes testes em ordem para confirmar que o ambiente está correto.

### 3.1 Testar o Docker

```bash
# Verificar a versão instalada
docker --version

# Rodar o container de teste oficial do Docker
docker run hello-world
```

Se o Docker estiver funcionando, você verá:
```
Hello from Docker!
This message shows that your installation appears to be working correctly.
```

O que aconteceu: o Docker baixou uma imagem mínima chamada `hello-world` do **Docker Hub** (repositório público de imagens em hub.docker.com) e a executou dentro de um container. Esse é exatamente o mecanismo que o Nextflow usará para cada ferramenta do pipeline.

### 3.2 Testar o Nextflow com Docker

```bash
# Rodar um pipeline de hello-world do Nextflow usando Docker
nextflow run hello -profile docker
```

Você deve ver:
```
N E X T F L O W  ~  version 24.x.x
executor >  local (4)
[xx/xxxxxx] process > sayHello (1) [100%] 4 of 4
Hola world!
Hello world!
Ciao world!
Bonjour world!
```

### 3.3 Testar o pipeline nf-core/rnaseq com dados de exemplo

Este é o teste mais importante: ele baixa e executa o pipeline completo com um pequeno dataset de exemplo.

```bash
# Criar diretório para o teste
mkdir -p ~/teste_rnaseq/resultados

# Rodar o pipeline com o perfil de teste
# AVISO: pode demorar 10-30 minutos na primeira vez
# (Nextflow baixa todos os containers Docker automaticamente)
nextflow run nf-core/rnaseq \
    -profile test,docker \
    --outdir ~/teste_rnaseq/resultados
```

**O que acontece durante essa execução:**
1. O Nextflow baixa o código do pipeline nf-core/rnaseq do GitHub
2. Para cada ferramenta (STAR, Salmon, FastQC, etc.), o Nextflow verifica se o container Docker já está no cache local
3. Se não estiver, o container é baixado automaticamente do Docker Hub
4. Cada ferramenta é executada dentro do seu próprio container isolado
5. Os containers são armazenados em cache -execuções futuras são muito mais rápidas

Se finalizar sem erros, você verá:
```
-[nf-core/rnaseq] Pipeline completed successfully -
```

**Isso confirma que Docker + Nextflow + nf-core/rnaseq estão funcionando corretamente.**

---

## Parte 4 - Baixando os dados do tutorial (SRA)

Vamos usar o dataset **Xia et al, 2021** (SRP282921), que é um experimento de QuantSeq 3' RNA-Seq em camundongos, sem UMIs. São 36 amostras disponíveis gratuitamente no NCBI SRA.

### 4.1 O que é o SRA e o nf-core/fetchngs

O **SRA (Sequence Read Archive)** é o maior repositório público de dados de sequenciamento do mundo, mantido pelo NCBI. Todos os artigos científicos publicados são obrigados a depositar os dados brutos lá.

O pipeline **nf-core/fetchngs** automatiza o download de dados do SRA. Como ele também usa Docker, não é necessário instalar nenhuma ferramenta de download manualmente -o Nextflow cuida de tudo.

O nf-core/fetchngs faz três coisas:
- Baixa os arquivos FASTQ diretamente do SRA
- Converte do formato SRA para FASTQ (quando necessário)
- Gera automaticamente um samplesheet compatível com nf-core/rnaseq

### 4.2 Preparar o diretório e o arquivo de IDs

O arquivo `ids.txt` com o ID do estudo já está neste repositório. Basta copiar para o seu diretório de trabalho:

```bash
# Criar o diretório do tutorial
mkdir -p ~/tutorial-quantseq
cd ~/tutorial-quantseq

# Criar o arquivo com o ID do estudo SRA
# (ou copiar o ids.txt que está neste repositório)
echo "SRP282921" > ids.txt

# Verificar o conteúdo
cat ids.txt
# Deve mostrar: SRP282921
```

**Por que um arquivo de IDs:** o nf-core/fetchngs aceita como entrada um arquivo de texto simples com os IDs do SRA que você quer baixar. Pode ser um ID de estudo (SRP...), de amostra (SRX...) ou de run (SRR...).

### 4.3 Baixar os dados com nf-core/fetchngs

```bash
mkdir -p ~/tutorial-quantseq/fastq

# Rodar o nf-core/fetchngs para baixar os FASTQs
# Tudo roda dentro de containers Docker -não precisa instalar nada
nextflow run nf-core/fetchngs \
    --input ids.txt \
    --outdir ~/tutorial-quantseq/fastq \
    -profile docker
```

O download de 36 amostras pode demorar bastante. Para rodar em segundo plano e salvar o log:

```bash
nohup nextflow run nf-core/fetchngs \
    --input ids.txt \
    --outdir ~/tutorial-quantseq/fastq \
    -profile docker \
    > ~/tutorial-quantseq/log_fetchngs.txt 2>&1 &

# Acompanhar o progresso em tempo real:
tail -f ~/tutorial-quantseq/log_fetchngs.txt
```

### 4.4 O que o fetchngs baixa

Após o download, você encontrará:

```
~/tutorial-quantseq/fastq/
├── fastq/
│   ├── SRR12345678_1.fastq.gz   # arquivo FASTQ de cada amostra
│   └── ...                       # (QuantSeq é single-end, só há R1)
└── samplesheet/
    └── samplesheet.csv           # samplesheet pronto para nf-core/rnaseq!
```

O `samplesheet.csv` gerado automaticamente já está no formato correto para ser usado diretamente na próxima etapa.

---

## Parte 5 - Preparando o samplesheet

O **samplesheet** é um arquivo CSV que informa ao pipeline: quais amostras processar, onde estão os arquivos FASTQ e qual é a orientação da biblioteca.

### 5.1 Samplesheet gerado automaticamente

Se você usou o nf-core/fetchngs, o samplesheet já está pronto em:
```
~/tutorial-quantseq/fastq/samplesheet/samplesheet.csv
```

Verifique o conteúdo:
```bash
head -5 ~/tutorial-quantseq/fastq/samplesheet/samplesheet.csv
```

### 5.2 Estrutura do samplesheet

O samplesheet tem 4 colunas obrigatórias:

| Coluna | O que é | Exemplo |
|--------|---------|---------|
| `sample` | Nome da amostra (único por amostra biológica) | `CTRL_REP1` |
| `fastq_1` | Caminho completo para o arquivo FASTQ R1 | `/home/user/fastq/SRR123_1.fastq.gz` |
| `fastq_2` | Caminho para FASTQ R2 (vazio se single-end) | _(vazio)_ |
| `strandedness` | Orientação da biblioteca | `forward` |

**Para QuantSeq:** os dados são **single-end** (só há R1, sem R2) e a biblioteca é **forward** (característica do protocolo QuantSeq da Lexogen).

### 5.3 Exemplo de samplesheet para QuantSeq

Se precisar criar ou ajustar manualmente (use o `samplesheet_template.csv` deste repositório como modelo):

```csv
sample,fastq_1,fastq_2,strandedness
CTRL_REP1,/home/usuario/tutorial-quantseq/fastq/SRR12345678_1.fastq.gz,,forward
CTRL_REP2,/home/usuario/tutorial-quantseq/fastq/SRR12345679_1.fastq.gz,,forward
CTRL_REP3,/home/usuario/tutorial-quantseq/fastq/SRR12345680_1.fastq.gz,,forward
TREAT_REP1,/home/usuario/tutorial-quantseq/fastq/SRR12345681_1.fastq.gz,,forward
TREAT_REP2,/home/usuario/tutorial-quantseq/fastq/SRR12345682_1.fastq.gz,,forward
TREAT_REP3,/home/usuario/tutorial-quantseq/fastq/SRR12345683_1.fastq.gz,,forward
```

**Regras importantes:**
- Use **caminhos absolutos** (começando com `/`), não relativos
- Para dados single-end, a coluna `fastq_2` fica **vazia** (mas a vírgula permanece)
- Não use espaços nos nomes das amostras (use underscores `_`)
- Os arquivos devem estar compactados com gzip (extensão `.fastq.gz`)

### 5.4 Baixar o genoma de referência

O dataset Xia et al usa camundongo. Vamos usar o genoma **GRCm38** (mm10) com a anotação **Gencode release M17**.

```bash
mkdir -p ~/tutorial-quantseq/referencia
cd ~/tutorial-quantseq/referencia

# Genoma do camundongo GRCm38 (FASTA)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M17/GRCm38.primary_assembly.genome.fa.gz

# Anotação do genoma (GTF)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M17/gencode.vM17.primary_assembly.annotation.gtf.gz

# Descompactar os dois arquivos
gunzip GRCm38.primary_assembly.genome.fa.gz
gunzip gencode.vM17.primary_assembly.annotation.gtf.gz

cd ~/tutorial-quantseq
```

**O que são esses dois arquivos:**
- **FASTA (`.fa`):** contém a sequência de DNA de todos os cromossomos do camundongo -é o "mapa" onde os reads serão alinhados
- **GTF (`.gtf`):** contém as coordenadas de todos os genes conhecidos -define onde cada gene começa e termina, necessário para quantificar a expressão

---

## Parte 6 - Executando o pipeline nf-core/rnaseq

### 6.1 Verificar se o Docker está ativo

```bash
sudo systemctl status docker
# Deve mostrar: Active: active (running)

# Se estiver parado:
sudo systemctl start docker
```

### 6.2 Configurar os parâmetros do pipeline

Crie o arquivo de parâmetros `params.yaml` (ou use o `params_exemplo.yaml` deste repositório como base):

```bash
cat > ~/tutorial-quantseq/params.yaml << 'EOF'
# Parâmetros do pipeline nf-core/rnaseq
# Dataset: Xia et al 2021 (SRP282921) - QuantSeq 3' RNA-Seq, camundongo

input: '/home/SEU_USUARIO/tutorial-quantseq/fastq/samplesheet/samplesheet.csv'
outdir: '/home/SEU_USUARIO/tutorial-quantseq/resultados_rnaseq'

fasta: '/home/SEU_USUARIO/tutorial-quantseq/referencia/GRCm38.primary_assembly.genome.fa'
gtf: '/home/SEU_USUARIO/tutorial-quantseq/referencia/gencode.vM17.primary_assembly.annotation.gtf'

aligner: 'star_salmon'

# Remove os primeiros 12 nucleotídeos (região de baixa qualidade no QuantSeq 3')
extra_trimgalore_args: '--clip_r1 12'
EOF
```

Substitua `SEU_USUARIO` pelo seu nome de usuário. Para saber qual é:
```bash
echo $USER
```

### 6.3 Executar o pipeline

```bash
cd ~/tutorial-quantseq

nextflow run nf-core/rnaseq \
    -profile docker \
    -params-file params.yaml \
    -r 3.10.1
```

**Explicação de cada argumento:**

| Argumento | Função |
|-----------|--------|
| `nf-core/rnaseq` | Nome do pipeline -o Nextflow baixa automaticamente do GitHub |
| `-profile docker` | Instrui o Nextflow a executar cada ferramenta dentro de um container Docker |
| `-params-file params.yaml` | Lê os parâmetros do arquivo que criamos |
| `-r 3.10.1` | Versão específica do pipeline -garante que o resultado é reprodutível |

### 6.4 O que acontece durante a execução

Você verá uma saída parecida com:
```
N E X T F L O W  ~  version 24.x.x
Launching `nf-core/rnaseq` [jolly_torvalds] DSL2 - revision: ...

executor >  local (142)
[xx/xxxxxx] process > NFCORE_RNASEQ:RNASEQ:PREPARE_GENOME:GTF_FILTER     [100%] 1 of 1
[xx/xxxxxx] process > NFCORE_RNASEQ:RNASEQ:PREPARE_GENOME:STAR_GENOMEGENERATE [  0%] 0 of 1
...
```

A primeira execução é mais lenta porque o Nextflow precisa baixar todos os containers Docker. Para um experimento com 36 amostras, espere entre 4 e 24 horas dependendo dos recursos do computador.

### 6.5 Rodar em segundo plano

Para não perder a execução caso o terminal feche:

```bash
# Opção 1: com nohup (mais simples)
nohup nextflow run nf-core/rnaseq \
    -profile docker \
    -params-file params.yaml \
    -r 3.10.1 \
    > log_rnaseq.txt 2>&1 &

# Acompanhar o progresso:
tail -f log_rnaseq.txt

# Opção 2: com screen (permite reconectar ao terminal)
screen -S rnaseq_tutorial

nextflow run nf-core/rnaseq \
    -profile docker \
    -params-file params.yaml \
    -r 3.10.1

# Para desanexar sem parar o pipeline: Ctrl+A, depois D
# Para reconectar: screen -r rnaseq_tutorial
```

### 6.6 Retomar uma execução interrompida

Se o pipeline falhar ou for interrompido (queda de energia, falta de memória, etc.), use `-resume`:

```bash
nextflow run nf-core/rnaseq \
    -profile docker \
    -params-file params.yaml \
    -r 3.10.1 \
    -resume
```

O `-resume` é uma das melhores funcionalidades do Nextflow: ele verifica quais etapas já foram concluídas e retoma de onde parou, sem refazer o que já foi feito. Isso economiza muito tempo.

---

## Parte 7 - Entendendo os resultados

Após a execução bem-sucedida, os resultados estarão em `~/tutorial-quantseq/resultados_rnaseq/`.

### 7.1 Estrutura dos resultados

```
resultados_rnaseq/
├── fastqc/                      # Qualidade dos reads brutos
│   └── *.html                   # Abrir no navegador
├── trimgalore/                  # Reads após limpeza dos adaptadores
├── star_salmon/                 # Alinhamento e quantificação
│   ├── salmon.merged.gene_counts.tsv       # <-- ESTE É O PRINCIPAL
│   ├── salmon.merged.gene_tpm.tsv          # Expressão em TPM
│   └── */                                  # Resultados por amostra
├── multiqc/
│   └── multiqc_report.html      # <-- COMECE AQUI para avaliar qualidade
└── pipeline_info/
    └── execution_report.html    # Relatório de execução do pipeline
```

### 7.2 Abrir o relatório MultiQC

O relatório **MultiQC** consolida todas as métricas de qualidade em um único arquivo HTML interativo. É o primeiro lugar a olhar depois da execução:

```bash
# Abrir no navegador (Ubuntu)
xdg-open resultados_rnaseq/multiqc/multiqc_report.html
```

**O que verificar no MultiQC:**

| Métrica | O que é | Valor aceitável |
|---------|---------|-----------------|
| FastQC Quality Scores | Qualidade dos reads (Phred score) | > 20 na maioria das posições |
| Trimming Stats | Reads retidos após limpeza | > 85% |
| STAR Alignment | Taxa de reads mapeados no genoma | > 70% para eucariotos |
| Gene Body Coverage | Distribuição de reads ao longo dos genes | Para QuantSeq: acumulação na extremidade 3' é esperado |
| PCA Plot | Agrupamento das amostras | Réplicas biológicas devem se agrupar |

### 7.3 A matriz de contagem de genes

O arquivo mais importante para análise downstream é:
```
resultados_rnaseq/star_salmon/salmon.merged.gene_counts.tsv
```

Este arquivo contém uma tabela onde:
- Cada **linha** é um gene
- Cada **coluna** é uma amostra
- Os **valores** são o número de reads que mapearam em cada gene (raw counts)

Esta matriz é a entrada para análises de expressão diferencial com DESeq2 ou edgeR.

```bash
# Ver as primeiras linhas da matriz de contagem
head -5 resultados_rnaseq/star_salmon/salmon.merged.gene_counts.tsv
```

### 7.4 Próximos passos

Com a matriz de contagem em mãos, os próximos passos típicos são:

1. **Análise de expressão diferencial (DEGs):** usando DESeq2 ou edgeR em R
2. **Análise de enriquecimento (GSEA/GO):** identificar quais funções biológicas estão enriquecidas
3. **Visualizações:** volcano plot, heatmap, PCA

Veja o README principal do projeto (`../README.md`) para instruções detalhadas de como fazer a análise com DESeq2.

---

## Solução de problemas comuns

### Docker não está rodando

```bash
# Verificar status
sudo systemctl status docker

# Iniciar
sudo systemctl start docker
```

### "permission denied while trying to connect to the Docker daemon"

Seu usuário não tem permissão para usar o Docker sem sudo.

```bash
sudo usermod -aG docker $USER
newgrp docker

# Testar
docker run hello-world
```

### "No space left on device"

O Docker acumulou imagens antigas que estão ocupando espaço em disco.

```bash
# Ver quanto espaço o Docker está usando
docker system df

# Remover imagens, containers e volumes não usados
docker system prune -a

# Verificar espaço em disco
df -h
```

### Pipeline muito lento ou sem resposta

```bash
# Ver os containers em execução
docker ps

# Ver uso de CPU e RAM em tempo real
docker stats
```

Se a RAM estiver esgotada, adicione ao `params.yaml`:
```yaml
max_memory: '12.GB'
max_cpus: 4
```

### "OutOfMemoryError" na etapa de indexação do STAR

O STAR precisa de ~30 GB de RAM para indexar o genoma do camundongo. Se não houver RAM suficiente, use o pseudo-alinhamento (usa muito menos memória):

```yaml
# Adicionar ao params.yaml:
pseudo_aligner: 'salmon'
skip_alignment: true
```

### Nextflow não encontra o samplesheet ou os FASTQs

Use sempre **caminhos absolutos**. Para descobrir o caminho absoluto de um arquivo:

```bash
realpath samplesheet.csv
```

### Pipeline falhou no meio, como retomar

```bash
nextflow run nf-core/rnaseq \
    -profile docker \
    -params-file params.yaml \
    -r 3.10.1 \
    -resume
```

### Ver os logs de um processo que falhou

O Nextflow informa o hash do processo que falhou (ex: `[ab/cd1234]`). Os logs ficam em:

```bash
# Log completo do Nextflow
cat .nextflow.log

# Logs do processo específico
cat work/ab/cd1234*/[nome_do_processo]/.command.log
cat work/ab/cd1234*/[nome_do_processo]/.command.err
```

---

## Dicionário rápido de termos

| Termo | Significado |
|-------|-------------|
| **FASTQ** | Formato de arquivo de sequenciamento (sequência + qualidade de cada base) |
| **Docker** | Plataforma de containers para executar software de forma isolada e reprodutível |
| **Container** | Ambiente isolado com um programa e todas as suas dependências |
| **Imagem Docker** | O "molde" a partir do qual containers são criados (armazenada no Docker Hub) |
| **Docker Hub** | Repositório público de imagens Docker (hub.docker.com) |
| **Nextflow** | Motor de workflow que orquestra a execução de pipelines bioinformáticos |
| **nf-core** | Comunidade que desenvolve e mantém pipelines Nextflow padronizados |
| **Pipeline** | Sequência automatizada de etapas de processamento bioinformático |
| **Samplesheet** | Arquivo CSV que lista as amostras e onde estão os FASTQs |
| **SRA** | Sequence Read Archive -repositório de dados de sequenciamento do NCBI |
| **GTF/GFF** | Formato de anotação genômica (coordenadas de genes, exons, etc.) |
| **FASTA** | Formato de arquivo de sequência de DNA/RNA |
| **STAR** | Alinhador de reads de RNA-Seq ao genoma |
| **Salmon** | Ferramenta de quantificação de expressão gênica |
| **MultiQC** | Ferramenta que consolida relatórios de qualidade em um HTML interativo |
| **TPM** | Transcripts Per Million -métrica de expressão normalizada |
| **Raw counts** | Contagem bruta de reads por gene (entrada para análise diferencial) |
| **DEGs** | Differentially Expressed Genes -genes diferencialmente expressos |

---

## Referências e leituras complementares

- **Tutorial original (em inglês):** https://tomsing1.github.io/blog/posts/nextflow-core-quantseq-1-settings/
- **Documentação nf-core/rnaseq:** https://nf-co.re/rnaseq
- **Documentação nf-core/fetchngs:** https://nf-co.re/fetchngs
- **Documentação do Nextflow:** https://www.nextflow.io/docs/latest/
- **Documentação do Docker:** https://docs.docker.com/engine/install/ubuntu/
- **Dataset Xia et al 2021 (SRP282921):** https://www.ncbi.nlm.nih.gov/sra/SRP282921
- **Dataset Nugent et al 2020 (SRP213880):** https://www.ncbi.nlm.nih.gov/sra/SRP213880
