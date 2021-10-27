{
 "cells": [
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# HSE21_HW1\n"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Создание символических ссылок на файлы\n",
    "ln -s /usr/share/data-minor-bioinf/assembly/oilMP_S4_L001_R1_001.fastq \n",
    "\n",
    "ln -s /usr/share/data-minor-bioinf/assembly/oilMP_S4_L001_R2_001.fastq \n",
    "\n",
    "ln -s /usr/share/data-minor-bioinf/assembly/oil_R1.fastq \n",
    "\n",
    "ln -s /usr/share/data-minor-bioinf/assembly/oil_R2.fastq "
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Выбираем образец 1024 и случайно отбираем 5000000 чтений pair-end и 1500000 чтений mate pair\n",
    "seqtk sample -s1024 oil_R1.fastq 5000000 > R1_paired_end.fastq\n",
    "\n",
    "seqtk sample -s1024 oil_R2.fastq 5000000 > R2_paired_end.fastq\n",
    "\n",
    "seqtk sample -s1024 oilMP_S4_L001_R1_001.fastq 1500000 > R1_mate_pairs.fastq\n",
    "\n",
    "seqtk sample -s1024 oilMP_S4_L001_R2_001.fastq 1500000 > R2_mate_pairs.fastq"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Оцениваем качество исходных чтений с помощью программ fastQC и multiQC, получаем общие данные \n",
    "mkdir fastq\n",
    "\n",
    "ls *.fastq | xargs -P 4 -tI{} fastqc -o fastqc {}\n",
    "\n",
    "mkdir multiqc\n",
    "\n",
    "multiqc -o multiqc fastqc"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Общая статистика для исходных чтений\n",
    "<img src=\"Desktop/до1.png\" width=\"800\"/>\n",
    "<img src=\"Desktop/до2.png\" width=\"800\"/>\n",
    "<img src=\"Desktop/до3.png\" width=\"800\"/>"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### С помощью программ platanus_trim и platanus_internal_trim подрезать чтения по качеству и удалить праймеры\n",
    "platanus_trim R1_paired_end.fastq R2_paired_end.fastq\n",
    "platanus_internal_trim R1_mate_pairs.fastq R2_mate_pairs.fastq"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Оцениваем качество \"подрезанных\" чтений c помощью программ fastQC и multiQС, получаем по ним общие данные\n",
    "mkdir trimmed_fastqc\n",
    "ls *trimmed | xargs -P 4 -tI{} fastqc -o trimmed_fastqc {}\n",
    "mkdir trimmed_multiqc\n",
    "multiqc -o trimmed_multiqc trimmed_fastqc"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Удаляем исходные файлы\n",
    "rm R1_paired_end.fastq\n",
    "\n",
    "rm R2_paired_end.fastq\n",
    "\n",
    "rm R1_mate_pairs.fastq\n",
    "\n",
    "rm R2_mate_pairs.fastq"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Общая статистика для подрезанных чтений\n",
    "<img src=\"Desktop/после1.png\" width=\"800\"/>\n",
    "<img src=\"Desktop/после2.png\" width=\"800\"/>\n",
    "<img src=\"Desktop/после3.png\" width=\"800\"/>"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Сравнив результаты для исходных и подрезанных чтений, можно сделать следующие выводы:\n",
    "1) Mean quality scores после удаления праймеров переместился в зеленую зону\n",
    "2) Уменьшились длины последовательностей\n",
    "3) Уменьшился процент содержания адаптеров"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Собираем контиги и скаффолды\n",
    "time platanus assemble -o Poil -t 2 -m 28 -f R1_paired_end.fastq.trimmed R2_paired_end.fastq.trimmed 2> assemble.log\n",
    "\n",
    "time platanus scaffold -o Poil -t 2 -c Poil_contig.fa -IP1 R1_paired_end.fastq.trimmed R2_paired_end.fastq.trimmed -OP2 R1_mate_pairs.fastq.int_trimmed R2_mate_pairs.fastq.int_trimmed 2> scaffold.log"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Код для анализа контигов"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "import numpy as np\n",
    "def analysis_c(contig):\n",
    "    counter = 0\n",
    "    length = 0\n",
    "    array = []\n",
    "    for row in contig:\n",
    "        if row[0] == '>':\n",
    "            counter += 1\n",
    "            length += int(row.split('_')[1][3:])\n",
    "            array.append(int(row.split('_')[1][3:]))\n",
    "    array = np.asarray(array)\n",
    "    array = np.sort(array)[::-1]\n",
    "    value = np.sum(array) / 2\n",
    "    n50 = array[np.cumsum(array) >= value][0]\n",
    "    \n",
    "    print('Общее количество контигов: ', counter)\n",
    "    print('Суммарная длина контигов: ', length)\n",
    "    print('Длина самого длинного контига: ', max(array))\n",
    "    print('N50: ', n50)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Код для анализа скаффолдов"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 3,
   "metadata": {},
   "outputs": [],
   "source": [
    "def analysis_s(scaffold):\n",
    "    counter = 0\n",
    "    length = 0\n",
    "    array = []\n",
    "    for row in scaffold:\n",
    "        if row[0] == '>':\n",
    "            counter += 1\n",
    "            length += int(row.split('_')[1][3:])\n",
    "            array.append(int(row.split('_')[1][3:]))\n",
    "    array = np.asarray(array)\n",
    "    array = np.sort(array)[::-1]\n",
    "    value = np.sum(array) / 2\n",
    "    n50 = array[np.cumsum(array) >= value][0]\n",
    "    \n",
    "    print('Общее количество скаффолдов: ', counter)\n",
    "    print('Суммарная длина скаффолдов: ', length)\n",
    "    print('Длина самого длинного скаффолда: ', max(array))\n",
    "    print('N50: ', n50)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Анализ контигов"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 7,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Общее количество контигов:  600\n",
      "Суммарная длина контигов:  3923404\n",
      "Длина самого длинного контига:  179304\n",
      "N50:  47798\n"
     ]
    }
   ],
   "source": [
    "contig_file = open('Poil_contig.fa', 'r')\n",
    "contig = contig_file.readlines()\n",
    "analysis_c(contig)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Анализ скаффолдов"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 8,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Общее количество скаффолдов:  70\n",
      "Суммарная длина скаффолдов:  3872710\n",
      "Длина самого длинного скаффолда:  3831215\n",
      "N50:  3831215\n"
     ]
    }
   ],
   "source": [
    "scaffold_file = open('Poil_scaffold.fa', 'r')\n",
    "scaffold = scaffold_file.readlines()\n",
    "analysis_s(scaffold)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Поиск наибольшего скаффолда"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 9,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      ">scaffold1_len3831215_cov232\n",
      "\n"
     ]
    }
   ],
   "source": [
    "for line in scaffold:\n",
    "    if int(line.split('_')[1][3:]) == 3831215:\n",
    "        print(line)\n",
    "        break"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.8.5"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 4
}
