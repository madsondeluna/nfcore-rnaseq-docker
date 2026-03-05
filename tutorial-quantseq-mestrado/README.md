# Tutorial Prático: RNA-Seq com nf-core/rnaseq e Docker

Tutorial completo para alunos de pós-graduação que estão aprendendo a montar e anotar dados de RNA-Seq do zero, utilizando o pipeline [nf-core/rnaseq](https://nf-co.re/rnaseq) dentro de containers Docker.

> Este tutorial é baseado no artigo: [QuantSeq RNAseq com nf-core (Thomas Singer, 2023)](https://tomsing1.github.io/blog/posts/nextflow-core-quantseq-1-settings/)

---

## O que você vai aprender

1. O que é Docker e por que ele é essencial para bioinformática reprodutível
2. Como instalar o Docker Desktop e usar o Docker pelo terminal
3. Como instalar o Nextflow e as ferramentas nf-core
4. Como verificar se tudo está funcionando corretamente
5. Como baixar os dados de sequenciamento públicos (FASTQs do SRA)
6. Como montar o samplesheet e configurar o pipeline
7. Como executar o pipeline nf-core/rnaseq com Docker
8. Como interpretar os resultados

---

## Sumário

1. [Contexto: o que é QuantSeq RNA-Seq](#1-contexto-o-que-é-quantseq-rna-seq)
2. [Por que usar Docker em bioinformática](#2-por-que-usar-docker-em-bioinformática)
3. [Parte 1 - Instalando o Docker](#parte-1---instalando-o-docker)
4. [Parte 2 - Instalando o Nextflow](#parte-2---instalando-o-nextflow)
5. [Parte 3 - Instalando as ferramentas nf-core](#parte-3---instalando-as-ferramentas-nf-core)
6. [Parte 4 - Verificando se tudo está funcionando](#parte-4---verificando-se-tudo-está-funcionando)
7. [Parte 5 - Baixando os dados do tutorial (SRA)](#parte-5---baixando-os-dados-do-tutorial-sra)
8. [Parte 6 - Preparando o samplesheet](#parte-6---preparando-o-samplesheet)
9. [Parte 7 - Executando o pipeline nf-core/rnaseq](#parte-7---executando-o-pipeline-nf-corernaseq)
10. [Parte 8 - Entendendo os resultados](#parte-8---entendendo-os-resultados)
11. [Solução de problemas comuns](#solução-de-problemas-comuns)

---

## 1. Contexto: o que é QuantSeq RNA-Seq

O **QuantSeq 3' mRNA-Seq** (da empresa Lexogen) é um protocolo de sequenciamento de RNA que captura apenas a região 3' dos transcritos. Isso reduz o custo e o tempo de processamento em comparação ao RNA-Seq convencional (que captura o transcrito inteiro), mantendo a capacidade de medir expressão gênica.

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

Imagine que você instala um programa de bioinformática hoje e ele funciona perfeitamente. Daqui a um ano, um colega tenta reproduzir sua análise no computador dele - mas o programa não funciona, porque a versão do Python é diferente, ou uma biblioteca foi atualizada, ou o sistema operacional é outro.

Esse problema de **reprodutibilidade** é um dos maiores desafios da bioinformática moderna.

### O que é um container Docker

O Docker resolve isso criando **containers**: pacotes completos e isolados que contêm o programa, todas as suas dependências, as versões exatas das bibliotecas e até o sistema operacional mínimo necessário para rodar. É como uma "caixa fechada" que funciona identicamente em qualquer computador que tenha Docker instalado.

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

O pipeline nf-core/rnaseq usa **dezenas de ferramentas diferentes** (STAR, Salmon, FastQC, Trim Galore, MultiQC, etc.). Em vez de você instalar cada uma manualmente (o que levaria horas e geraria conflitos de versão), o Nextflow baixa automaticamente os containers Docker de cada ferramenta quando necessário.

Você só precisa ter o Docker instalado. O resto é automático.

---

## Parte 1 - Instalando o Docker

### Opção A: Docker Desktop (recomendado para iniciantes)

O **Docker Desktop** é uma aplicação gráfica com interface visual que instala tudo automaticamente. É a forma mais fácil de começar.

#### macOS (Apple Silicon M1/M2/M3 ou Intel)

1. Acesse: https://www.docker.com/products/docker-desktop/
2. Clique em **"Download Docker Desktop"**
3. Escolha a versão correta:
   - **Apple Silicon** se seu Mac tem chip M1, M2 ou M3
   - **Intel** se seu Mac tem processador Intel (Macs mais antigos)
4. Abra o arquivo `.dmg` baixado
5. Arraste o ícone do Docker para a pasta **Applications**
6. Abra o Docker Desktop pelo Launchpad ou pela pasta Applications
7. Aguarde a inicialização (um ícone de baleia aparecerá na barra de menu no topo da tela)
8. Quando a baleia parar de se mover, o Docker está pronto

**Como verificar se o Mac é Apple Silicon ou Intel:**
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
4. **IMPORTANTE:** durante a instalação, aceite a opção de habilitar o **WSL 2** (Windows Subsystem for Linux)
   - O WSL 2 é necessário para que o Docker funcione bem no Windows
   - Se solicitado, siga as instruções para instalar o WSL 2 pelo Windows Update
5. Reinicie o computador se solicitado
6. Abra o Docker Desktop pelo menu Iniciar
7. Aguarde a inicialização completa

#### Linux (Ubuntu/Debian)

No Linux, o Docker Desktop também está disponível, mas muitos usuários preferem instalar apenas o motor Docker (Docker Engine) pelo terminal:

```bash
# 1. Remover versões antigas (se houver)
sudo apt-get remove docker docker-engine docker.io containerd runc

# 2. Atualizar os pacotes do sistema
sudo apt-get update

# 3. Instalar dependências necessárias
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

# 5. Adicionar o repositório oficial do Docker
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.gpg] \
  https://download.docker.com/linux/ubuntu \
  $(lsb_release -cs) stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null

# 6. Instalar o Docker Engine
sudo apt-get update
sudo apt-get install -y docker-ce docker-ce-cli containerd.io docker-compose-plugin
```

**Por que cada passo:** os passos 1-3 preparam o sistema, o passo 4 adiciona a chave de segurança do Docker (evita instalar software modificado/malicioso), o passo 5 registra o repositório oficial, e o passo 6 instala o Docker Engine de fato.

**Configuração pós-instalação no Linux** (importante para não precisar de `sudo` toda hora):

```bash
# Adicionar seu usuário ao grupo docker
# Isso permite rodar comandos docker sem sudo
sudo usermod -aG docker $USER

# Aplicar a mudança sem precisar fazer logout
newgrp docker

# Iniciar o serviço Docker automaticamente ao ligar o computador
sudo systemctl enable docker
sudo systemctl start docker
```

---

### Opção B: Docker via terminal (linha de comando pura)

Se você usa macOS e prefere instalar pelo terminal, pode usar o **Homebrew**:

```bash
# Instalar o Homebrew primeiro (se não tiver)
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Instalar o Docker Desktop via Homebrew Cask
brew install --cask docker

# Abrir o Docker Desktop para inicializar
open /Applications/Docker.app
```

---

### Como ativar/iniciar o Docker

O Docker precisa estar **em execução** para o Nextflow conseguir criar containers. Veja como verificar e iniciar:

**macOS e Windows:**
- Abra o aplicativo **Docker Desktop** pela pasta Applications (macOS) ou menu Iniciar (Windows)
- Aguarde o ícone da baleia na barra de menu ficar estável (para de se mover)
- O Docker está pronto quando você vê "Docker Desktop is running" no painel do aplicativo

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

O **Nextflow** é o motor de workflow que orquestra a execução do pipeline. Ele é quem lê as instruções do nf-core/rnaseq, baixa os containers Docker de cada ferramenta e executa tudo na ordem correta.

### 2.1 Verificar se o Java está instalado

O Nextflow precisa do Java para funcionar. Verifique:

```bash
java -version
```

**O que você deve ver (algo parecido com):**
```
openjdk version "17.0.9" 2023-10-17
OpenJDK Runtime Environment Temurin-17.0.9+9 (build 17.0.9+9)
OpenJDK Server VM Temurin-17.0.9+9 (build 17.0.9+9, mixed mode)
```

Se o Java não estiver instalado ou for anterior à versão 11, instale:

```bash
# Ubuntu/Debian
sudo apt-get install -y default-jdk

# macOS com Homebrew
brew install openjdk@17

# Após instalar no macOS, adicionar ao PATH:
echo 'export PATH="/opt/homebrew/opt/openjdk@17/bin:$PATH"' >> ~/.zshrc
source ~/.zshrc
```

### 2.2 Instalar o Nextflow

```bash
# Baixar o instalador do Nextflow
curl -s https://get.nextflow.io | bash

# O comando acima cria um arquivo executável chamado "nextflow"
# no diretório atual. Mova para um local acessível em todo o sistema:

# macOS/Linux:
chmod +x nextflow
sudo mv nextflow /usr/local/bin/

# Verificar a instalação
nextflow -version
```

**O que você deve ver:**
```
      N E X T F L O W
      version 24.x.x build ...
      created ...
      cite doi:10.1038/nbt.3820
      http://nextflow.io
```

**Por que mover para `/usr/local/bin/`:** esse diretório é automaticamente incluído no PATH do sistema, o que permite você chamar `nextflow` de qualquer pasta, sem precisar escrever o caminho completo.

---

## Parte 3 - Instalando as ferramentas nf-core

O **nf-core** é uma comunidade que mantém pipelines bioinformáticos padronizados. Além dos pipelines em si (que o Nextflow baixa automaticamente), existe um pacote de ferramentas de linha de comando que facilita tarefas como baixar pipelines, gerar samplesheets e baixar dados do SRA.

### 3.1 Instalar as ferramentas nf-core via pip

```bash
# Verificar se o pip está instalado
pip --version

# Se não estiver, instalar:
# Ubuntu/Debian:
sudo apt-get install -y python3-pip

# macOS:
brew install python3

# Instalar as ferramentas nf-core
pip install nf-core

# Ou, se preferir instalar apenas para o seu usuário (sem sudo):
pip install --user nf-core
```

### 3.2 Verificar a instalação do nf-core

```bash
nf-core --version
```

**O que você deve ver:**
```
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 2.x.x - https://nf-co.re
```

---

## Parte 4 - Verificando se tudo está funcionando

Antes de rodar seus dados reais, execute estes testes para confirmar que o ambiente está correto.

### 4.1 Verificar Docker

```bash
# Verificar a versão do Docker
docker --version
# Esperado: Docker version 24.x.x, build ...

# Testar se o Docker consegue baixar e rodar containers
docker run hello-world
```

**O que você deve ver após `docker run hello-world`:**
```
Hello from Docker!
This message shows that your installation appears to be working correctly.
...
```

Se aparecer essa mensagem, o Docker está funcionando. O que aconteceu: o Docker baixou uma imagem mínima chamada `hello-world` do Docker Hub (repositório público de imagens) e a executou dentro de um container.

### 4.2 Verificar Nextflow com Docker

```bash
# Testar o Nextflow com um pipeline simples de exemplo
# Este comando baixa e executa um hello-world do Nextflow
nextflow run hello -profile docker
```

**O que deve acontecer:** o Nextflow vai baixar um pipeline de teste, criar containers Docker para cada processo e executar. No final, você deve ver algo como:
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

O nf-core/rnaseq inclui um perfil de teste (`-profile test`) que usa dados minúsculos para verificar se o pipeline funciona:

```bash
# Criar um diretório para o teste
mkdir -p ~/teste_rnaseq/resultados

# Rodar o pipeline com o perfil de teste
# Isso vai baixar o pipeline e todos os containers automaticamente
# AVISO: pode demorar 10-30 minutos na primeira vez (baixando containers)
nextflow run nf-core/rnaseq \
    -profile test,docker \
    --outdir ~/teste_rnaseq/resultados
```

**O que acontece durante a execução:**
- O Nextflow baixa o código do pipeline nf-core/rnaseq do GitHub
- Para cada ferramenta (STAR, Salmon, FastQC, etc.), o Nextflow baixa automaticamente o container Docker correspondente do Docker Hub
- Os containers são armazenados localmente em cache, então da segunda vez é muito mais rápido
- O pipeline executa com um pequeno dataset de teste

**Se finalizar sem erros**, você verá:
```
-[nf-core/rnaseq] Pipeline completed successfully -
```

---

## Parte 5 - Baixando os dados do tutorial (SRA)

Vamos usar o dataset **Xia et al, 2021** (SRP282921), que é um experimento de QuantSeq 3' RNA-Seq em camundongos, sem UMIs. São 36 amostras disponíveis gratuitamente no NCBI SRA.

### 5.1 O que é o SRA e o nf-core/fetchngs

O **SRA (Sequence Read Archive)** é o maior repositório público de dados de sequenciamento do mundo, mantido pelo NCBI. Todos os artigos publicados são obrigados a depositar os dados brutos lá.

O pipeline **nf-core/fetchngs** automatiza o download de dados do SRA, incluindo:
- Download dos arquivos FASTQ
- Geração automática de um samplesheet compatível com nf-core/rnaseq
- Conversão do formato SRA para FASTQ

### 5.2 Criar o arquivo de IDs para download

```bash
# Criar o diretório do tutorial
mkdir -p ~/tutorial-quantseq
cd ~/tutorial-quantseq

# Criar o arquivo com o ID do estudo SRA
echo "SRP282921" > ids.txt

# Verificar o conteúdo
cat ids.txt
```

**Por que um arquivo de IDs:** o nf-core/fetchngs aceita como entrada um arquivo de texto simples com os IDs do SRA que você quer baixar. Pode ser um ID de estudo (SRP...), de amostra (SRX...) ou de run (SRR...).

### 5.3 Baixar os dados com nf-core/fetchngs

```bash
# Criar os diretórios
mkdir -p ~/tutorial-quantseq/{fastq,resultados_rnaseq}

# Rodar o nf-core/fetchngs para baixar os FASTQs
nextflow run nf-core/fetchngs \
    --input ids.txt \
    --outdir ~/tutorial-quantseq/fastq \
    -profile docker

# O download pode demorar bastante dependendo da sua conexão.
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

Após o download, você encontrará:
```
~/tutorial-quantseq/fastq/
├── fastq/
│   ├── SRR12345678_1.fastq.gz   # arquivo R1 (forward)
│   ├── SRR12345678_2.fastq.gz   # arquivo R2 (reverse) -- se paired-end
│   └── ...                       # (QuantSeq é single-end, só terá R1)
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

# Baixar uma amostra específica (exemplo)
fastq-dump --split-files --gzip SRR12345678

# Para baixar todas as amostras de um estudo, use o nf-core/fetchngs (mais fácil)
```

---

## Parte 6 - Preparando o samplesheet

O **samplesheet** é um arquivo CSV que informa ao pipeline: quais amostras processar, onde estão os arquivos FASTQ, e qual é a orientação da biblioteca.

### 6.1 Samplesheet gerado automaticamente pelo fetchngs

Se você usou o nf-core/fetchngs, o samplesheet já está pronto em:
```
~/tutorial-quantseq/fastq/samplesheet/samplesheet.csv
```

Verifique o conteúdo:
```bash
head -5 ~/tutorial-quantseq/fastq/samplesheet/samplesheet.csv
```

### 6.2 Estrutura do samplesheet

O samplesheet tem 4 colunas obrigatórias:

```csv
sample,fastq_1,fastq_2,strandedness
```

| Coluna | O que é | Exemplo |
|--------|---------|---------|
| `sample` | Nome da amostra (único por amostra biológica) | `CTRL_REP1` |
| `fastq_1` | Caminho completo para o arquivo FASTQ R1 | `/home/user/fastq/SRR123_1.fastq.gz` |
| `fastq_2` | Caminho para FASTQ R2 (vazio se single-end) | _(vazio)_ |
| `strandedness` | Orientação da biblioteca | `forward` |

**Para QuantSeq:** os dados são **single-end** (só há R1) e a biblioteca é **forward** (strand-specific sentido forward, característica do protocolo QuantSeq da Lexogen).

### 6.3 Exemplo de samplesheet para QuantSeq

Se precisar criar manualmente, o formato é:

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
- Use **caminhos absolutos** (começando com `/`), não relativos
- Para dados single-end, a coluna `fastq_2` fica **vazia** (mas a vírgula permanece)
- Não use espaços nos nomes das amostras (use underscores `_`)
- Os arquivos devem estar compactados com gzip (extensão `.fastq.gz`)

### 6.4 Verificar o genoma de referência

O dataset Xia et al usa camundongo. Vamos usar o genoma **GRCm38** (mm10) com a anotação **Gencode release M17**.

Para baixar o genoma de referência:

```bash
mkdir -p ~/tutorial-quantseq/referencia
cd ~/tutorial-quantseq/referencia

# Genoma do camundongo GRCm38 (FASTA)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M17/GRCm38.primary_assembly.genome.fa.gz

# Anotação do genoma (GTF)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M17/gencode.vM17.primary_assembly.annotation.gtf.gz

# Descompactar
gunzip GRCm38.primary_assembly.genome.fa.gz
gunzip gencode.vM17.primary_assembly.annotation.gtf.gz
```

**Por que esses arquivos:**
- **FASTA (`.fa`):** contém a sequência de DNA de todos os cromossomos do camundongo - é o "mapa" onde os reads serão alinhados
- **GTF (`.gtf`):** contém as coordenadas de todos os genes conhecidos - define onde cada gene começa e termina, necessário para quantificar a expressão

---

## Parte 7 - Executando o pipeline nf-core/rnaseq

### 7.1 Configurar os parâmetros do pipeline

Crie o arquivo de parâmetros `params.yaml`:

```bash
cat > ~/tutorial-quantseq/params.yaml << 'EOF'
# Parâmetros do pipeline nf-core/rnaseq
# Dataset: Xia et al 2021 (SRP282921) - QuantSeq 3' RNA-Seq, camundongo

# Arquivos de entrada
input: '/home/SEU_USUARIO/tutorial-quantseq/fastq/samplesheet/samplesheet.csv'
outdir: '/home/SEU_USUARIO/tutorial-quantseq/resultados_rnaseq'

# Genoma de referência
fasta: '/home/SEU_USUARIO/tutorial-quantseq/referencia/GRCm38.primary_assembly.genome.fa'
gtf: '/home/SEU_USUARIO/tutorial-quantseq/referencia/gencode.vM17.primary_assembly.annotation.gtf'

# Alinhador: STAR + Salmon (padrão e recomendado)
aligner: 'star_salmon'

# Configurações específicas para QuantSeq 3' mRNA-Seq
# QuantSeq captura apenas a região 3', então precisamos
# recortar os primeiros 12 bases (poli-A da extremidade 3')
extra_trimgalore_args: '--clip_r1 12'

# Limitar memória e CPU se necessário (ajuste conforme seu computador)
# max_memory: '16.GB'
# max_cpus: 8
EOF
```

**Substituir `SEU_USUARIO`** pelo seu nome de usuário no sistema. Para saber qual é:
```bash
echo $USER
# ou
whoami
```

### 7.2 Verificar se o Docker está rodando

```bash
# Verificar se o Docker está ativo
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

**Explicação de cada parâmetro:**
- `nf-core/rnaseq`: nome do pipeline (o Nextflow baixa automaticamente do GitHub)
- `-profile docker`: instrui o Nextflow a usar Docker para cada ferramenta
- `-params-file params.yaml`: lê os parâmetros do arquivo que criamos
- `-r 3.10.1`: versão específica do pipeline (garante reprodutibilidade)

### 7.4 O que acontece durante a execução

A primeira execução é mais lenta porque o Nextflow precisa baixar todos os containers Docker. Para um experimento com 36 amostras, espere entre 4 e 24 horas dependendo do computador.

Você verá uma saída parecida com:
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

**Retomando uma execução interrompida** (queda de energia, timeout, etc.):

```bash
nextflow run nf-core/rnaseq \
    -profile docker \
    -params-file params.yaml \
    -r 3.10.1 \
    -resume
```

O `-resume` é uma das melhores funcionalidades do Nextflow: ele verifica quais etapas já foram concluídas e retoma de onde parou, sem refazer o que já foi feito.

### 7.5 Rodar em segundo plano

Para execuções longas, use `screen` ou `nohup` para que o pipeline continue mesmo se você fechar o terminal:

```bash
# Opção 1: com nohup (mais simples)
nohup nextflow run nf-core/rnaseq \
    -profile docker \
    -params-file params.yaml \
    -r 3.10.1 \
    > log_rnaseq.txt 2>&1 &

# Ver o progresso:
tail -f log_rnaseq.txt

# Opção 2: com screen (mais controle)
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

Após a execução bem-sucedida, os resultados estarão em `~/tutorial-quantseq/resultados_rnaseq/`.

### 8.1 Estrutura dos resultados

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

### 8.2 Primeiro: abrir o MultiQC

O relatório **MultiQC** consolida todas as métricas de qualidade em um único arquivo HTML interativo. É o primeiro lugar a olhar depois da execução:

```bash
# macOS
open resultados_rnaseq/multiqc/multiqc_report.html

# Linux
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

### 8.3 A matriz de contagem de genes

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

### 8.4 Próximos passos

Com a matriz de contagem em mãos, os próximos passos típicos são:

1. **Análise de expressão diferencial (DEGs):** usando DESeq2 ou edgeR em R
2. **Análise de enriquecimento (GSEA/GO):** identificar quais funções biológicas estão enriquecidas
3. **Visualizações:** volcano plot, heatmap, PCA

Veja o README principal do projeto (`../README.md`) para instruções detalhadas de como fazer a análise com DESeq2.

---

## Solução de problemas comuns

### "Docker daemon is not running" ou "Cannot connect to Docker daemon"

O Docker não está iniciado.

```bash
# macOS/Windows: abra o aplicativo Docker Desktop e aguarde inicializar

# Linux:
sudo systemctl start docker
sudo systemctl status docker   # deve mostrar "active (running)"
```

### "permission denied while trying to connect to the Docker daemon"

Seu usuário não tem permissão para usar o Docker sem sudo.

```bash
# Adicionar usuário ao grupo docker
sudo usermod -aG docker $USER

# Aplicar sem precisar de logout (ou faça logout/login)
newgrp docker

# Testar
docker run hello-world
```

### "No space left on device"

O Docker acumulou muitas imagens antigas. Libere espaço:

```bash
# Ver quanto espaço o Docker está usando
docker system df

# Remover imagens, containers e volumes não usados
docker system prune -a

# Verificar espaço em disco
df -h
```

### Pipeline muito lento ou travado

```bash
# Ver os processos em execução no Docker
docker ps

# Ver uso de recursos
docker stats

# Se necessário, limitar recursos no params.yaml:
# max_memory: '8.GB'
# max_cpus: 4
```

### "OutOfMemoryError" ou processo morto por falta de RAM

A etapa de indexação do STAR precisa de muita RAM (~30 GB para o genoma humano/camundongo). Se seu computador tem menos RAM:

```bash
# Adicionar ao params.yaml:
# max_memory: '12.GB'

# Ou usar o pseudo-alinhamento (Salmon sem STAR) que usa muito menos RAM:
# pseudo_aligner: 'salmon'
# skip_alignment: true
```

### Nextflow não encontra o arquivo de samplesheet

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

# Os diretórios de trabalho ficam em ./work/
# O Nextflow informa o hash do processo que falhou, ex: [ab/cd1234]
# Os logs desse processo ficam em:
cat work/ab/cd1234.../[nome_do_processo]/.command.log
cat work/ab/cd1234.../[nome_do_processo]/.command.err
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
| **SRA** | Sequence Read Archive - repositório de dados de sequenciamento do NCBI |
| **GTF/GFF** | Formato de anotação genômica (coordenadas de genes, exons, etc.) |
| **FASTA** | Formato de arquivo de sequência de DNA/RNA |
| **STAR** | Alinhador de reads de RNA-Seq ao genoma |
| **Salmon** | Ferramenta de quantificação de expressão gênica |
| **MultiQC** | Ferramenta que consolida relatórios de qualidade em um HTML interativo |
| **TPM** | Transcripts Per Million - métrica de expressão normalizada |
| **Raw counts** | Contagem bruta de reads por gene (entrada para análise diferencial) |
| **DEGs** | Differentially Expressed Genes - genes diferencialmente expressos |

---

## Referências e leituras complementares

- **Tutorial original (em inglês):** https://tomsing1.github.io/blog/posts/nextflow-core-quantseq-1-settings/
- **Documentação nf-core/rnaseq:** https://nf-co.re/rnaseq
- **Documentação nf-core/fetchngs:** https://nf-co.re/fetchngs
- **Documentação do Nextflow:** https://www.nextflow.io/docs/latest/
- **Documentação do Docker:** https://docs.docker.com/get-started/
- **Dataset Xia et al 2021 (SRP282921):** https://www.ncbi.nlm.nih.gov/sra/SRP282921
- **Dataset Nugent et al 2020 (SRP213880):** https://www.ncbi.nlm.nih.gov/sra/SRP213880
