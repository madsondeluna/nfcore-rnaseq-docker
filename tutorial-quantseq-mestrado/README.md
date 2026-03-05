# Tutorial Pratico: RNA-Seq com nf-core/rnaseq e Docker

Tutorial completo para alunos de pos-graduacao que estao aprendendo a montar e anotar dados de RNA-Seq do zero, utilizando o pipeline [nf-core/rnaseq](https://nf-co.re/rnaseq) dentro de containers Docker.

> Este tutorial e baseado no artigo: [QuantSeq RNAseq com nf-core (Thomas Singer, 2023)](https://tomsing1.github.io/blog/posts/nextflow-core-quantseq-1-settings/)

---

## O que voce vai aprender

1. O que e Docker e por que ele e essencial para bioinformatica reprodutivel
2. Como instalar o Docker Desktop e usar o Docker pelo terminal
3. Como instalar o Nextflow e as ferramentas nf-core
4. Como verificar se tudo esta funcionando corretamente
5. Como baixar os dados de sequenciamento publicos (FASTQs do SRA)
6. Como montar o samplesheet e configurar o pipeline
7. Como executar o pipeline nf-core/rnaseq com Docker
8. Como interpretar os resultados

---

## Sumario

1. [Contexto: o que e QuantSeq RNA-Seq](#1-contexto-o-que-e-quantseq-rna-seq)
2. [Por que usar Docker em bioinformatica](#2-por-que-usar-docker-em-bioinformatica)
3. [Parte 1 - Instalando o Docker](#parte-1---instalando-o-docker)
4. [Parte 2 - Instalando o Nextflow](#parte-2---instalando-o-nextflow)
5. [Parte 3 - Instalando as ferramentas nf-core](#parte-3---instalando-as-ferramentas-nf-core)
6. [Parte 4 - Verificando se tudo esta funcionando](#parte-4---verificando-se-tudo-esta-funcionando)
7. [Parte 5 - Baixando os dados do tutorial (SRA)](#parte-5---baixando-os-dados-do-tutorial-sra)
8. [Parte 6 - Preparando o samplesheet](#parte-6---preparando-o-samplesheet)
9. [Parte 7 - Executando o pipeline nf-core/rnaseq](#parte-7---executando-o-pipeline-nf-corernaseq)
10. [Parte 8 - Entendendo os resultados](#parte-8---entendendo-os-resultados)
11. [Solucao de problemas comuns](#solucao-de-problemas-comuns)

---

## 1. Contexto: o que e QuantSeq RNA-Seq

O **QuantSeq 3' mRNA-Seq** (da empresa Lexogen) e um protocolo de sequenciamento de RNA que captura apenas a regiao 3' dos transcritos. Isso reduz o custo e o tempo de processamento em comparacao ao RNA-Seq convencional (que captura o transcrito inteiro), mantendo a capacidade de medir expressao genica.

**Caracteristicas importantes:**
- Dados **single-end** (apenas uma direcao de leitura, sem R2)
- Biblioteca **strand-specific** (sentido forward)
- Reads mais curtos (50-100 bp) concentrados na extremidade 3' dos genes
- Pode conter **UMIs** (Unique Molecular Identifiers) em alguns protocolos

**Dados que vamos usar neste tutorial:**

O tutorial original usa dois conjuntos de dados publicos de camundongo disponveis no NCBI SRA:

| Dataset | Artigo | Acesso SRA | Acesso GEO | Detalhes |
|---------|--------|-----------|-----------|----------|
| Xia et al, 2021 | Sem UMIs | SRP282921 | GSE158152 | 36 amostras, camundongo |
| Nugent et al, 2020 | Com UMIs | SRP213880 | GSE134031 | Camundongo |

Nos vamos baixar e processar o dataset **Xia et al (SRP282921)**, que e o mais simples (sem UMIs).

---

## 2. Por que usar Docker em bioinformatica

### O problema que o Docker resolve

Imagine que voce instala um programa de bioinformatica hoje e ele funciona perfeitamente. Daqui a um ano, um colega tenta reproduzir sua analise no computador dele - mas o programa nao funciona, porque a versao do Python e diferente, ou uma biblioteca foi atualizada, ou o sistema operacional e outro.

Esse problema de **reproducibilidade** e um dos maiores desafios da bioinformatica moderna.

### O que e um container Docker

O Docker resolve isso criando **containers**: pacotes completos e isolados que contem o programa, todas as suas dependencias, as versoes exatas das bibliotecas e ate o sistema operacional minimo necessario para rodar. E como uma "caixa fechada" que funciona identicamente em qualquer computador que tenha Docker instalado.

```
Sem Docker:                  Com Docker:
+-----------------+          +---------------------------+
| Seu computador  |          | Container Docker          |
| Python 3.11     |  ====>   | Python 3.9                |
| NumPy 1.24      |          | NumPy 1.21                |
| STAR 2.7.10a    |          | STAR 2.7.9a               |
| ...             |          | TUDO que o programa precisa|
+-----------------+          +---------------------------+
                             | Funciona igual em qualquer|
                             | computador com Docker!    |
                             +---------------------------+
```

### Por que o nf-core usa Docker

O pipeline nf-core/rnaseq usa **dezenas de ferramentas diferentes** (STAR, Salmon, FastQC, Trim Galore, MultiQC, etc.). Em vez de voce instalar cada uma manualmente (o que levaria horas e geraria conflitos de versao), o Nextflow baixa automaticamente os containers Docker de cada ferramenta quando necessario.

Voce so precisa ter o Docker instalado. O resto e automatico.

---

## Parte 1 - Instalando o Docker

### Opcao A: Docker Desktop (recomendado para iniciantes)

O **Docker Desktop** e uma aplicacao grafica com interface visual que instala tudo automaticamente. E a forma mais facil de comecar.

#### macOS (Apple Silicon M1/M2/M3 ou Intel)

1. Acesse: https://www.docker.com/products/docker-desktop/
2. Clique em **"Download Docker Desktop"**
3. Escolha a versao correta:
   - **Apple Silicon** se seu Mac tem chip M1, M2 ou M3
   - **Intel** se seu Mac tem processador Intel (Macs mais antigos)
4. Abra o arquivo `.dmg` baixado
5. Arraste o icone do Docker para a pasta **Applications**
6. Abra o Docker Desktop pelo Launchpad ou pela pasta Applications
7. Aguarde a inicializacao (um icone de baleia aparecera na barra de menu no topo da tela)
8. Quando a baleia parar de se mover, o Docker esta pronto

**Como verificar se o Mac e Apple Silicon ou Intel:**
```bash
# Abra o Terminal e execute:
uname -m

# Se aparecer "arm64" => Apple Silicon (M1/M2/M3)
# Se aparecer "x86_64" => Intel
```

#### Windows 10/11

1. Acesse: https://www.docker.com/products/docker-desktop/
2. Baixe o instalador para Windows
3. Execute o instalador `.exe`
4. **IMPORTANTE:** durante a instalacao, aceite a opcao de habilitar o **WSL 2** (Windows Subsystem for Linux)
   - O WSL 2 e necessario para que o Docker funcione bem no Windows
   - Se solicitado, siga as instrucoes para instalar o WSL 2 pelo Windows Update
5. Reinicie o computador se solicitado
6. Abra o Docker Desktop pelo menu Iniciar
7. Aguarde a inicializacao completa

#### Linux (Ubuntu/Debian)

No Linux, o Docker Desktop tambem esta disponivel, mas muitos usuarios preferem instalar apenas o motor Docker (Docker Engine) pelo terminal:

```bash
# 1. Remover versoes antigas (se houver)
sudo apt-get remove docker docker-engine docker.io containerd runc

# 2. Atualizar os pacotes do sistema
sudo apt-get update

# 3. Instalar dependencias necessarias
sudo apt-get install -y \
    ca-certificates \
    curl \
    gnupg \
    lsb-release

# 4. Adicionar a chave GPG oficial do Docker
# (garante que o software vem da fonte oficial)
sudo mkdir -p /etc/apt/keyrings
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | \
    sudo gpg --dearmor -o /etc/apt/keyrings/docker.gpg

# 5. Adicionar o repositorio oficial do Docker
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.gpg] \
  https://download.docker.com/linux/ubuntu \
  $(lsb_release -cs) stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null

# 6. Instalar o Docker Engine
sudo apt-get update
sudo apt-get install -y docker-ce docker-ce-cli containerd.io docker-compose-plugin
```

**Por que cada passo:** os passos 1-3 preparam o sistema, o passo 4 adiciona a chave de segurança do Docker (evita instalar software modificado/malicioso), o passo 5 registra o repositorio oficial, e o passo 6 instala o Docker Engine de fato.

**Configuracao pos-instalacao no Linux** (importante para nao precisar de `sudo` toda hora):

```bash
# Adicionar seu usuario ao grupo docker
# Isso permite rodar comandos docker sem sudo
sudo usermod -aG docker $USER

# Aplicar a mudanca sem precisar fazer logout
newgrp docker

# Iniciar o servico Docker automaticamente ao ligar o computador
sudo systemctl enable docker
sudo systemctl start docker
```

---

### Opcao B: Docker via terminal (linha de comando pura)

Se voce usa macOS e prefere instalar pelo terminal, pode usar o **Homebrew**:

```bash
# Instalar o Homebrew primeiro (se nao tiver)
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Instalar o Docker Desktop via Homebrew Cask
brew install --cask docker

# Abrir o Docker Desktop para inicializar
open /Applications/Docker.app
```

---

### Como ativar/iniciar o Docker

O Docker precisa estar **em execucao** para o Nextflow conseguir criar containers. Veja como verificar e iniciar:

**macOS e Windows:**
- Abra o aplicativo **Docker Desktop** pela pasta Applications (macOS) ou menu Iniciar (Windows)
- Aguarde o icone da baleia na barra de menu ficar estavel (para de se mover)
- O Docker esta pronto quando voce ve "Docker Desktop is running" no painel do aplicativo

**Linux:**
```bash
# Verificar o status do Docker
sudo systemctl status docker

# Se estiver parado, iniciar
sudo systemctl start docker

# Para iniciar automaticamente ao ligar o computador
sudo systemctl enable docker
```

---

## Parte 2 - Instalando o Nextflow

O **Nextflow** e o motor de workflow que orquestra a execucao do pipeline. Ele e quem le as instrucoes do nf-core/rnaseq, baixa os containers Docker de cada ferramenta e executa tudo na ordem correta.

### 2.1 Verificar se o Java esta instalado

O Nextflow precisa do Java para funcionar. Verifique:

```bash
java -version
```

**O que voce deve ver (algo parecido com):**
```
openjdk version "17.0.9" 2023-10-17
OpenJDK Runtime Environment Temurin-17.0.9+9 (build 17.0.9+9)
OpenJDK Server VM Temurin-17.0.9+9 (build 17.0.9+9, mixed mode)
```

Se o Java nao estiver instalado ou for anterior a versao 11, instale:

```bash
# Ubuntu/Debian
sudo apt-get install -y default-jdk

# macOS com Homebrew
brew install openjdk@17

# Apos instalar no macOS, adicionar ao PATH:
echo 'export PATH="/opt/homebrew/opt/openjdk@17/bin:$PATH"' >> ~/.zshrc
source ~/.zshrc
```

### 2.2 Instalar o Nextflow

```bash
# Baixar o instalador do Nextflow
curl -s https://get.nextflow.io | bash

# O comando acima cria um arquivo executavel chamado "nextflow"
# no diretorio atual. Mova para um local acessivel em todo o sistema:

# macOS/Linux:
chmod +x nextflow
sudo mv nextflow /usr/local/bin/

# Verificar a instalacao
nextflow -version
```

**O que voce deve ver:**
```
      N E X T F L O W
      version 24.x.x build ...
      created ...
      cite doi:10.1038/nbt.3820
      http://nextflow.io
```

**Por que mover para `/usr/local/bin/`:** esse diretorio e automaticamente incluido no PATH do sistema, o que permite voce chamar `nextflow` de qualquer pasta, sem precisar escrever o caminho completo.

---

## Parte 3 - Instalando as ferramentas nf-core

O **nf-core** e uma comunidade que mantem pipelines bioinformaticos padronizados. Alem dos pipelines em si (que o Nextflow baixa automaticamente), existe um pacote de ferramentas de linha de comando que facilita tarefas como baixar pipelines, gerar samplesheets e baixar dados do SRA.

### 3.1 Instalar as ferramentas nf-core via pip

```bash
# Verificar se o pip esta instalado
pip --version

# Se nao estiver, instalar:
# Ubuntu/Debian:
sudo apt-get install -y python3-pip

# macOS:
brew install python3

# Instalar as ferramentas nf-core
pip install nf-core

# Ou, se preferir instalar apenas para o seu usuario (sem sudo):
pip install --user nf-core
```

### 3.2 Verificar a instalacao do nf-core

```bash
nf-core --version
```

**O que voce deve ver:**
```
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 2.x.x - https://nf-co.re
```

---

## Parte 4 - Verificando se tudo esta funcionando

Antes de rodar seus dados reais, execute estes testes para confirmar que o ambiente esta correto.

### 4.1 Verificar Docker

```bash
# Verificar a versao do Docker
docker --version
# Esperado: Docker version 24.x.x, build ...

# Testar se o Docker consegue baixar e rodar containers
docker run hello-world
```

**O que voce deve ver apos `docker run hello-world`:**
```
Hello from Docker!
This message shows that your installation appears to be working correctly.
...
```

Se aparecer essa mensagem, o Docker esta funcionando. O que aconteceu: o Docker baixou uma imagem minima chamada `hello-world` do Docker Hub (repositorio publico de imagens) e a executou dentro de um container.

### 4.2 Verificar Nextflow com Docker

```bash
# Testar o Nextflow com um pipeline simples de exemplo
# Este comando baixa e executa um hello-world do Nextflow
nextflow run hello -profile docker
```

**O que deve acontecer:** o Nextflow vai baixar um pipeline de teste, criar containers Docker para cada processo e executar. No final, voce deve ver algo como:
```
N E X T F L O W  ~  version 24.x.x
Launching `https://github.com/nextflow-io/hello` [happy_...] DSL2 - revision: ...
executor >  local (4)
[xx/xxxxxx] process > sayHello (1) [100%] 4 of 4
Hola world!
Hello world!
Ciao world!
Bonjour world!
```

### 4.3 Testar o pipeline nf-core/rnaseq com dados de exemplo

O nf-core/rnaseq inclui um perfil de teste (`-profile test`) que usa dados minusculos para verificar se o pipeline funciona:

```bash
# Criar um diretorio para o teste
mkdir -p ~/teste_rnaseq/resultados

# Rodar o pipeline com o perfil de teste
# Isso vai baixar o pipeline e todos os containers automaticamente
# AVISO: pode demorar 10-30 minutos na primeira vez (baixando containers)
nextflow run nf-core/rnaseq \
    -profile test,docker \
    --outdir ~/teste_rnaseq/resultados
```

**O que acontece durante a execucao:**
- O Nextflow baixa o codigo do pipeline nf-core/rnaseq do GitHub
- Para cada ferramenta (STAR, Salmon, FastQC, etc.), o Nextflow baixa automaticamente o container Docker correspondente do Docker Hub
- Os containers sao armazenados localmente em cache, entao da segunda vez e muito mais rapido
- O pipeline executa com um pequeno dataset de teste

**Se finalizar sem erros**, voce vera:
```
-[nf-core/rnaseq] Pipeline completed successfully -
```

---

## Parte 5 - Baixando os dados do tutorial (SRA)

Vamos usar o dataset **Xia et al, 2021** (SRP282921), que e um experimento de QuantSeq 3' RNA-Seq em camundongos, sem UMIs. Sao 36 amostras disponveis gratuitamente no NCBI SRA.

### 5.1 O que e o SRA e o nf-core/fetchngs

O **SRA (Sequence Read Archive)** e o maior repositorio publico de dados de sequenciamento do mundo, mantido pelo NCBI. Todos os artigos publicados sao obrigados a depositar os dados brutos la.

O pipeline **nf-core/fetchngs** automatiza o download de dados do SRA, incluindo:
- Download dos arquivos FASTQ
- Geracao automatica de um samplesheet compativel com nf-core/rnaseq
- Conversao do formato SRA para FASTQ

### 5.2 Criar o arquivo de IDs para download

```bash
# Criar o diretorio do tutorial
mkdir -p ~/tutorial-quantseq
cd ~/tutorial-quantseq

# Criar o arquivo com o ID do estudo SRA
echo "SRP282921" > ids.txt

# Verificar o conteudo
cat ids.txt
```

**Por que um arquivo de IDs:** o nf-core/fetchngs aceita como entrada um arquivo de texto simples com os IDs do SRA que voce quer baixar. Pode ser um ID de estudo (SRP...), de amostra (SRX...) ou de run (SRR...).

### 5.3 Baixar os dados com nf-core/fetchngs

```bash
# Criar os diretorios
mkdir -p ~/tutorial-quantseq/{fastq,resultados_rnaseq}

# Rodar o nf-core/fetchngs para baixar os FASTQs
nextflow run nf-core/fetchngs \
    --input ids.txt \
    --outdir ~/tutorial-quantseq/fastq \
    -profile docker

# O download pode demorar bastante dependendo da sua conexao.
# Para rodar em segundo plano e salvar o log:
nohup nextflow run nf-core/fetchngs \
    --input ids.txt \
    --outdir ~/tutorial-quantseq/fastq \
    -profile docker \
    > ~/tutorial-quantseq/log_fetchngs.txt 2>&1 &

# Acompanhar o progresso:
tail -f ~/tutorial-quantseq/log_fetchngs.txt
```

**O que o nf-core/fetchngs baixa:**

Apos o download, voce encontrara:
```
~/tutorial-quantseq/fastq/
├── fastq/
│   ├── SRR12345678_1.fastq.gz   # arquivo R1 (forward)
│   ├── SRR12345678_2.fastq.gz   # arquivo R2 (reverse) -- se paired-end
│   └── ...                       # (QuantSeq e single-end, so tera R1)
└── samplesheet/
    └── samplesheet.csv           # samplesheet pronto para nf-core/rnaseq!
```

**Dica importante:** o nf-core/fetchngs gera automaticamente o `samplesheet.csv` no formato correto para ser usado diretamente com o nf-core/rnaseq. Isso economiza muito trabalho!

### 5.4 (Alternativa) Baixar manualmente via SRA-toolkit

Se preferir baixar manualmente, instale o SRA-toolkit:

```bash
# Ubuntu/Debian
sudo apt-get install -y sra-toolkit

# macOS com Homebrew
brew install sratoolkit

# Baixar uma amostra especifica (exemplo)
fastq-dump --split-files --gzip SRR12345678

# Para baixar todas as amostras de um estudo, use o nf-core/fetchngs (mais facil)
```

---

## Parte 6 - Preparando o samplesheet

O **samplesheet** e um arquivo CSV que informa ao pipeline: quais amostras processar, onde estao os arquivos FASTQ, e qual e a orientacao da biblioteca.

### 6.1 Samplesheet gerado automaticamente pelo fetchngs

Se voce usou o nf-core/fetchngs, o samplesheet ja esta pronto em:
```
~/tutorial-quantseq/fastq/samplesheet/samplesheet.csv
```

Verifique o conteudo:
```bash
head -5 ~/tutorial-quantseq/fastq/samplesheet/samplesheet.csv
```

### 6.2 Estrutura do samplesheet

O samplesheet tem 4 colunas obrigatorias:

```csv
sample,fastq_1,fastq_2,strandedness
```

| Coluna | O que e | Exemplo |
|--------|---------|---------|
| `sample` | Nome da amostra (unico por amostra biologica) | `CTRL_REP1` |
| `fastq_1` | Caminho completo para o arquivo FASTQ R1 | `/home/user/fastq/SRR123_1.fastq.gz` |
| `fastq_2` | Caminho para FASTQ R2 (vazio se single-end) | _(vazio)_ |
| `strandedness` | Orientacao da biblioteca | `forward` |

**Para QuantSeq:** os dados sao **single-end** (so ha R1) e a biblioteca e **forward** (strand-specific sentido forward, caracteristica do protocolo QuantSeq da Lexogen).

### 6.3 Exemplo de samplesheet para QuantSeq

Se precisar criar manualmente, o formato e:

```csv
sample,fastq_1,fastq_2,strandedness
CTRL_REP1,/home/user/tutorial-quantseq/fastq/SRR12345678_1.fastq.gz,,forward
CTRL_REP2,/home/user/tutorial-quantseq/fastq/SRR12345679_1.fastq.gz,,forward
CTRL_REP3,/home/user/tutorial-quantseq/fastq/SRR12345680_1.fastq.gz,,forward
TREAT_REP1,/home/user/tutorial-quantseq/fastq/SRR12345681_1.fastq.gz,,forward
TREAT_REP2,/home/user/tutorial-quantseq/fastq/SRR12345682_1.fastq.gz,,forward
TREAT_REP3,/home/user/tutorial-quantseq/fastq/SRR12345683_1.fastq.gz,,forward
```

**Regras importantes:**
- Use **caminhos absolutos** (comecando com `/`), nao relativos
- Para dados single-end, a coluna `fastq_2` fica **vazia** (mas a virgula permanece)
- Nao use espacos nos nomes das amostras (use underscores `_`)
- Os arquivos devem estar compactados com gzip (extensao `.fastq.gz`)

### 6.4 Verificar o genoma de referencia

O dataset Xia et al usa camundongo. Vamos usar o genoma **GRCm38** (mm10) com a anotacao **Gencode release M17**.

Para baixar o genoma de referencia:

```bash
mkdir -p ~/tutorial-quantseq/referencia
cd ~/tutorial-quantseq/referencia

# Genoma do camundongo GRCm38 (FASTA)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M17/GRCm38.primary_assembly.genome.fa.gz

# Anotacao do genoma (GTF)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M17/gencode.vM17.primary_assembly.annotation.gtf.gz

# Descompactar
gunzip GRCm38.primary_assembly.genome.fa.gz
gunzip gencode.vM17.primary_assembly.annotation.gtf.gz
```

**Por que esses arquivos:**
- **FASTA (`.fa`):** contem a sequencia de DNA de todos os cromossomos do camundongo - e o "mapa" onde os reads serao alinhados
- **GTF (`.gtf`):** contem as coordenadas de todos os genes conhecidos - define onde cada gene comeca e termina, necessario para quantificar a expressao

---

## Parte 7 - Executando o pipeline nf-core/rnaseq

### 7.1 Configurar os parametros do pipeline

Crie o arquivo de parametros `params.yaml`:

```bash
cat > ~/tutorial-quantseq/params.yaml << 'EOF'
# Parametros do pipeline nf-core/rnaseq
# Dataset: Xia et al 2021 (SRP282921) - QuantSeq 3' RNA-Seq, camundongo

# Arquivos de entrada
input: '/home/SEU_USUARIO/tutorial-quantseq/fastq/samplesheet/samplesheet.csv'
outdir: '/home/SEU_USUARIO/tutorial-quantseq/resultados_rnaseq'

# Genoma de referencia
fasta: '/home/SEU_USUARIO/tutorial-quantseq/referencia/GRCm38.primary_assembly.genome.fa'
gtf: '/home/SEU_USUARIO/tutorial-quantseq/referencia/gencode.vM17.primary_assembly.annotation.gtf'

# Alinhador: STAR + Salmon (padrao e recomendado)
aligner: 'star_salmon'

# Configuracoes especificas para QuantSeq 3' mRNA-Seq
# QuantSeq captura apenas a regiao 3', entao precisamos
# recortar os primeiros 12 bases (poli-A da extremidade 3')
extra_trimgalore_args: '--clip_r1 12'

# Limitar memoria e CPU se necessario (ajuste conforme seu computador)
# max_memory: '16.GB'
# max_cpus: 8
EOF
```

**Substituir `SEU_USUARIO`** pelo seu nome de usuario no sistema. Para saber qual e:
```bash
echo $USER
# ou
whoami
```

### 7.2 Verificar se o Docker esta rodando

```bash
# Verificar se o Docker esta ativo
docker info

# Se aparecer erro "Cannot connect to the Docker daemon":
# - macOS/Windows: abra o aplicativo Docker Desktop
# - Linux: sudo systemctl start docker
```

### 7.3 Rodar o pipeline

```bash
cd ~/tutorial-quantseq

# Executar o pipeline
nextflow run nf-core/rnaseq \
    -profile docker \
    -params-file params.yaml \
    -r 3.10.1
```

**Explicacao de cada parametro:**
- `nf-core/rnaseq`: nome do pipeline (o Nextflow baixa automaticamente do GitHub)
- `-profile docker`: instrui o Nextflow a usar Docker para cada ferramenta
- `-params-file params.yaml`: le os parametros do arquivo que criamos
- `-r 3.10.1`: versao especifica do pipeline (garante reproducibilidade)

### 7.4 O que acontece durante a execucao

A primeira execucao e mais lenta porque o Nextflow precisa baixar todos os containers Docker. Para um experimento com 36 amostras, espere entre 4 e 24 horas dependendo do computador.

Voce vera uma saida parecida com:
```
N E X T F L O W  ~  version 24.x.x
Launching `nf-core/rnaseq` [jolly_torvalds] DSL2 - revision: ...
Core Nextflow options
  runName        : jolly_torvalds
  containerEngine: docker
  ...

executor >  local (142)
[xx/xxxxxx] process > NFCORE_RNASEQ:RNASEQ:PREPARE_GENOME:GTF_FILTER     [100%] 1 of 1
[xx/xxxxxx] process > NFCORE_RNASEQ:RNASEQ:PREPARE_GENOME:STAR_GENOMEGENERATE [  0%] 0 of 1
...
```

**Retomando uma execucao interrompida** (queda de energia, timeout, etc.):

```bash
nextflow run nf-core/rnaseq \
    -profile docker \
    -params-file params.yaml \
    -r 3.10.1 \
    -resume
```

O `-resume` e uma das melhores funcionalidades do Nextflow: ele verifica quais etapas ja foram concluidas e retoma de onde parou, sem refazer o que ja foi feito.

### 7.5 Rodar em segundo plano

Para execucoes longas, use `screen` ou `nohup` para que o pipeline continue mesmo se voce fechar o terminal:

```bash
# Opcao 1: com nohup (mais simples)
nohup nextflow run nf-core/rnaseq \
    -profile docker \
    -params-file params.yaml \
    -r 3.10.1 \
    > log_rnaseq.txt 2>&1 &

# Ver o progresso:
tail -f log_rnaseq.txt

# Opcao 2: com screen (mais controle)
screen -S rnaseq_tutorial

nextflow run nf-core/rnaseq \
    -profile docker \
    -params-file params.yaml \
    -r 3.10.1

# Para desanexar (pipeline continua rodando): pressione Ctrl+A, depois D
# Para reconectar: screen -r rnaseq_tutorial
```

---

## Parte 8 - Entendendo os resultados

Apos a execucao bem-sucedida, os resultados estarao em `~/tutorial-quantseq/resultados_rnaseq/`.

### 8.1 Estrutura dos resultados

```
resultados_rnaseq/
├── fastqc/                      # Qualidade dos reads brutos
│   └── *.html                   # Abrir no navegador
├── trimgalore/                  # Reads apos limpeza dos adaptadores
├── star_salmon/                 # Alinhamento e quantificacao
│   ├── salmon.merged.gene_counts.tsv       # <-- ESTE E O PRINCIPAL
│   ├── salmon.merged.gene_tpm.tsv          # Expressao em TPM
│   └── */                                  # Resultados por amostra
├── multiqc/
│   └── multiqc_report.html      # <-- COMECE AQUI para avaliar qualidade
└── pipeline_info/
    └── execution_report.html    # Relatorio de execucao do pipeline
```

### 8.2 Primeiro: abrir o MultiQC

O relatorio **MultiQC** consolida todas as metricas de qualidade em um unico arquivo HTML interativo. E o primeiro lugar a olhar depois da execucao:

```bash
# macOS
open resultados_rnaseq/multiqc/multiqc_report.html

# Linux
xdg-open resultados_rnaseq/multiqc/multiqc_report.html
```

**O que verificar no MultiQC:**

| Metrica | O que e | Valor aceitavel |
|---------|---------|-----------------|
| FastQC Quality Scores | Qualidade dos reads (Phred score) | > 20 na maioria das posicoes |
| Trimming Stats | Reads retidos apos limpeza | > 85% |
| STAR Alignment | Taxa de reads mapeados no genoma | > 70% para eucariotos |
| Gene Body Coverage | Distribuicao de reads ao longo dos genes | Para QuantSeq: acumulacao na extremidade 3' e esperado |
| PCA Plot | Agrupamento das amostras | Replicas biologicas devem se agrupar |

### 8.3 A matriz de contagem de genes

O arquivo mais importante para analise downstream e:
```
resultados_rnaseq/star_salmon/salmon.merged.gene_counts.tsv
```

Este arquivo contem uma tabela onde:
- Cada **linha** e um gene
- Cada **coluna** e uma amostra
- Os **valores** sao o numero de reads que mapearam em cada gene (raw counts)

Esta matriz e a entrada para analises de expressao diferencial com DESeq2 ou edgeR.

```bash
# Ver as primeiras linhas da matriz de contagem
head -5 resultados_rnaseq/star_salmon/salmon.merged.gene_counts.tsv
```

### 8.4 Proximos passos

Com a matriz de contagem em maos, os proximos passos tipicos sao:

1. **Analise de expressao diferencial (DEGs):** usando DESeq2 ou edgeR em R
2. **Analise de enriquecimento (GSEA/GO):** identificar quais funcoes biologicas estao enriquecidas
3. **Visualizacoes:** volcano plot, heatmap, PCA

Veja o README principal do projeto (`../README.md`) para instrucoes detalhadas de como fazer a analise com DESeq2.

---

## Solucao de problemas comuns

### "Docker daemon is not running" ou "Cannot connect to Docker daemon"

O Docker nao esta iniciado.

```bash
# macOS/Windows: abra o aplicativo Docker Desktop e aguarde inicializar

# Linux:
sudo systemctl start docker
sudo systemctl status docker   # deve mostrar "active (running)"
```

### "permission denied while trying to connect to the Docker daemon"

Seu usuario nao tem permissao para usar o Docker sem sudo.

```bash
# Adicionar usuario ao grupo docker
sudo usermod -aG docker $USER

# Aplicar sem precisar de logout (ou faca logout/login)
newgrp docker

# Testar
docker run hello-world
```

### "No space left on device"

O Docker acumulou muitas imagens antigas. Libere espaco:

```bash
# Ver quanto espaco o Docker esta usando
docker system df

# Remover imagens, containers e volumes nao usados
docker system prune -a

# Verificar espaco em disco
df -h
```

### Pipeline muito lento ou travado

```bash
# Ver os processos em execucao no Docker
docker ps

# Ver uso de recursos
docker stats

# Se necessario, limitar recursos no params.yaml:
# max_memory: '8.GB'
# max_cpus: 4
```

### "OutOfMemoryError" ou processo morto por falta de RAM

A etapa de indexacao do STAR precisa de muita RAM (~30 GB para o genoma humano/camundongo). Se seu computador tem menos RAM:

```bash
# Adicionar ao params.yaml:
# max_memory: '12.GB'

# Ou usar o pseudo-alinhamento (Salmon sem STAR) que usa muito menos RAM:
# pseudo_aligner: 'salmon'
# skip_alignment: true
```

### Nextflow nao encontra o arquivo de samplesheet

Use sempre **caminhos absolutos** no samplesheet e no params.yaml:

```bash
# Para descobrir o caminho absoluto de um arquivo:
realpath samplesheet.csv
```

### Pipeline falhou no meio, como retomar

```bash
# SEMPRE use -resume para retomar
nextflow run nf-core/rnaseq \
    -profile docker \
    -params-file params.yaml \
    -r 3.10.1 \
    -resume
```

### Ver os logs detalhados de um processo que falhou

```bash
# Ver os logs do Nextflow
cat .nextflow.log

# Os diretorios de trabalho ficam em ./work/
# O Nextflow informa o hash do processo que falhou, ex: [ab/cd1234]
# Os logs desse processo ficam em:
cat work/ab/cd1234.../[nome_do_processo]/.command.log
cat work/ab/cd1234.../[nome_do_processo]/.command.err
```

---

## Dicionario rapido de termos

| Termo | Significado |
|-------|-------------|
| **FASTQ** | Formato de arquivo de sequenciamento (sequencia + qualidade de cada base) |
| **Docker** | Plataforma de containers para executar software de forma isolada e reprodutivel |
| **Container** | Ambiente isolado com um programa e todas as suas dependencias |
| **Imagem Docker** | O "molde" a partir do qual containers sao criados (armazenada no Docker Hub) |
| **Docker Hub** | Repositorio publico de imagens Docker (hub.docker.com) |
| **Nextflow** | Motor de workflow que orquestra a execucao de pipelines bioinformaticos |
| **nf-core** | Comunidade que desenvolve e mantem pipelines Nextflow padronizados |
| **Pipeline** | Sequencia automatizada de etapas de processamento bioinformatico |
| **Samplesheet** | Arquivo CSV que lista as amostras e onde estao os FASTQs |
| **SRA** | Sequence Read Archive - repositorio de dados de sequenciamento do NCBI |
| **GTF/GFF** | Formato de anotacao genomica (coordenadas de genes, exons, etc.) |
| **FASTA** | Formato de arquivo de sequencia de DNA/RNA |
| **STAR** | Alinhador de reads de RNA-Seq ao genoma |
| **Salmon** | Ferramenta de quantificacao de expressao genica |
| **MultiQC** | Ferramenta que consolida relatorios de qualidade em um HTML interativo |
| **TPM** | Transcripts Per Million - metrica de expressao normalizada |
| **Raw counts** | Contagem bruta de reads por gene (entrada para analise diferencial) |
| **DEGs** | Differentially Expressed Genes - genes diferencialmente expressos |

---

## Referencias e leituras complementares

- **Tutorial original (em ingles):** https://tomsing1.github.io/blog/posts/nextflow-core-quantseq-1-settings/
- **Documentacao nf-core/rnaseq:** https://nf-co.re/rnaseq
- **Documentacao nf-core/fetchngs:** https://nf-co.re/fetchngs
- **Documentacao do Nextflow:** https://www.nextflow.io/docs/latest/
- **Documentacao do Docker:** https://docs.docker.com/get-started/
- **Dataset Xia et al 2021 (SRP282921):** https://www.ncbi.nlm.nih.gov/sra/SRP282921
- **Dataset Nugent et al 2020 (SRP213880):** https://www.ncbi.nlm.nih.gov/sra/SRP213880
