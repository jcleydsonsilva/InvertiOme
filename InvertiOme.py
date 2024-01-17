#!/usr/bin/env python
# coding: utf-8

# In[1]:


import subprocess
import os
import sys
from Bio import SeqIO

# jupyter nbconvert --to script nome_do_notebook.ipynb

def run_local_blast(query_sequence, subject_sequence,window_size,perc_identity,gqry,gsubj,outfmt=6):
    # Criar um arquivo temporário para armazenar as sequências
    query_file = "query.fasta"
    subject_file = "subject.fasta"
    output_file = query_sequence[0].replace('>','').replace('-','')
    output_file += '_' + subject_sequence[0].replace('>','').replace('-','') + "_results.txt"
    
    
    gquery_file = "gquery.fasta"
    gsubject_file = "gsubject.fasta"
    goutput_file = query_sequence[0].replace('>','').replace('-','')
    goutput_file += 'g_' + subject_sequence[0].replace('>','').replace('-','') + "g_results.txt"    

    
    output_pdf = query_sequence[0].replace('>','').replace('-','')
    output_pdf += '_' + subject_sequence[0].replace('>','').replace('-','') + "_results.pdf"

    goutput_pdf = query_sequence[0].replace('>','').replace('-','')
    goutput_pdf += '_' + subject_sequence[0].replace('>','').replace('-','') + "_Gene_results.pdf"
    
    with open(query_file, "w") as f:
        f.write(">"+query_sequence[0]+"\n" + query_sequence[1]+"\n")
    f.close()
    
    with open(subject_file, "w") as f:
        f.write(">"+subject_sequence[0]+"\n" + subject_sequence[1]+"\n")
    f.close()
    
    with open(gquery_file, "w") as f:
        for key in gqry:
            f.write(f">{key}\n{gqry[key]}\n")
    f.close()

    with open(gsubject_file, "w") as f:
        for key in gsubj:
            f.write(f">{key}\n{gsubj[key]}\n")
    f.close()
    
    # Executar o BLAST localmente
    command = f"makeblastdb -in {subject_file} -dbtype nucl \n"
    command += f"blastn -num_threads 4 -dust yes -word_size {window_size} -perc_identity {perc_identity} -evalue 0.0000000001 -query {query_file} -db {subject_file} -outfmt {outfmt} -out {output_file}\n"   
    command += f"makeblastdb -in {gsubject_file} -dbtype nucl \n"
    
    command += f"blastn -qcov_hsp_perc 60 -num_alignments 10 -dust yes -num_threads 4 -perc_identity 0.6 -evalue 0.0000000000001 -query {gquery_file} -db {gsubject_file} -outfmt {outfmt} -out {goutput_file}\n"
    command += f"cat  {goutput_file} | sed 's/_/\t/g' > {goutput_file}_f.txt\n"
    
    scpt = open('script.R','w')
    scpt.write(f"pdf('{output_pdf}')\n")
    scpt.write(f"blastnData = read.table ('{output_file}', sep = '\t', header=FALSE)\n")
    scpt.write(f"resolucao_MB <- 10000000\n")
    scpt.write(f"escala_x <- seq(0, max(blastnData$V9), by = resolucao_MB)\n")
    scpt.write(f"escala_y <- seq(0, max(blastnData$V7), by = resolucao_MB)\n")
    scpt.write(f"escala_x_rotulos <- escala_x / 1000000\n")
    scpt.write(f"escala_y_rotulos <- escala_y / 1000000\n")
    scpt.write(f"options(scipen = 999)\n")
    scpt.write(f"cores <- ifelse(blastnData$V10 - blastnData$V9 < 0, 'red', 'black')\n")
    scpt.write(f"plot (blastnData$V9, blastnData$V7, cex = .250, col = cores,xaxt = 'n', yaxt = 'n', main = '', xlab = '{query_sequence[0]}', ylab = '{subject_sequence[0]}')\n")
    scpt.write(f"axis(1, at = escala_x, labels = escala_x_rotulos, las = 2)\n")
    scpt.write(f"axis(2, at = escala_y, labels = escala_y_rotulos, las = 2)\n")
    scpt.write(f"dev.off()\n")
    scpt.close()

    scpt = open('script_Gene.R','w')
    scpt.write(f"pdf('{goutput_pdf}')\n")
    scpt.write(f"blastnData = read.table ('{output_file}', sep = '\t', header=FALSE)\n")
    scpt.write(f"gblastnData = read.table ('{goutput_file}_f.txt', sep = '\t', header=FALSE)\n")
    scpt.write(f"resolucao_MB <- 10000000\n")
    scpt.write(f"escala_x <- seq(0, max(gblastnData$V4), by = resolucao_MB)\n")
    scpt.write(f"escala_y <- seq(0, max(gblastnData$V2), by = resolucao_MB)\n")
    scpt.write(f"escala_x_rotulos <- escala_x / 1000000\n")
    scpt.write(f"escala_y_rotulos <- escala_y / 1000000\n")
    scpt.write(f"options(scipen = 999)\n")
    scpt.write(f"gcores <- ifelse(gblastnData$V12 - gblastnData$V11 < 0, 'green', 'blue')\n")
    scpt.write(f"plot (gblastnData$V4, gblastnData$V2, cex = .250, col = gcores,xaxt = 'n', yaxt = 'n', main = '', xlab = '{query_sequence[0]}', ylab = '{subject_sequence[0]}')\n")
    scpt.write(f"axis(1, at = escala_x, labels = escala_x_rotulos, las = 2)\n")
    scpt.write(f"axis(2, at = escala_y, labels = escala_y_rotulos, las = 2)\n")
    scpt.write(f"dev.off()\n")
    scpt.close()
    
    command += f"Rscript script.R\n"
    command += f"Rscript script_Gene.R\n"
    cmd = open('cmd.bash','w')
    cmd.write(command)
    cmd.close()
    
    cmd = f"bash cmd.bash"
    
    print (command)
    process = subprocess.Popen(cmd,stdout=subprocess.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()
    # Ler os resultados do BLAST
    with open(output_file, "r") as f:
        blast_results = f.read()

    # Remover os arquivos temporários
    os.remove(query_file)
    os.remove(subject_file)
    os.remove(output_file)
    return blast_results

def plot_blast_hits(blast_results):
    # Seu código de plotagem aqui
    pass

def main():
    import json
    name_f = "config.json"
    with open(name_f, 'r') as files:
        json = json.load(files)
    files.close()
    
    #files = list()
    #f = open(json['fastafiles'],'r')
    #for l in f:
    #    files.append(l.strip().split('\t'))
    #f.close()
    
    gff3 = {}
    f = open(json['gff3files'],'r')
    for l in f:
        print (l.split('-')[1])
        gff3[l.split('-')[1]] = l.strip()
    f.close()
    
    # Ler a sequência do cromossomo 4 do arquivo FASTA
    chrs = {}
    for record in SeqIO.parse(json['fastafiles'], "fasta"):
        chrs[record.id] = record.seq
    
    i = 0
    for f1 in chrs:
        j = 0
        for f2 in chrs:            
            if f1 != f2 and i < j:
                kf1 = f1.replace('>','').split('_')[0] # B73v5_chr10
                kf2 = f2.replace('>','').split('_')[0] # B73v5_chr10
                #
                gqry = {}
                genesqry = open(gff3[kf1],'r')
                for l in genesqry:
                    if l.startswith(str(json['chr']).lower()):
                        l = l.split()
                        if l[2].startswith('gene'):
                            key = f"{l[8].split(';')[0].replace('ID=','')}_{l[3]}"
                            gqry[key] = chrs[f1][int(l[3]):int(l[4])]
                genesqry.close()
                #
                gsubj = {}
                genessubj = open(gff3[kf2],'r')
                for l in genessubj:
                    if l.startswith(str(json['chr']).lower()):
                        l = l.split()
                        if l[2].startswith('gene'):
                            key = f"{l[8].split(';')[0].replace('ID=','')}_{l[3]}"
                            gsubj[key] = chrs[f2][int(l[3]):int(l[4])]
                genessubj.close()
                
                #for i in range(int(json['start']), int(json['end']), int(json['window_size'])):
                # Definir o tamanho da janela
                window_size = int(json['word_size'])

                # Executar o BLAST para cada janela
                qry = list()
                qry.append(f1)
                qry.append(str(chrs[f1]))
                #qry.append(chr4_qry[int(i):i+int(json['window_size']) ])

                subj = list()
                subj.append(f2)
                subj.append(str(chrs[f2]))
                #subj.append( chr4_subj[int(json['start']):int(json['end'])])

                blast_results = run_local_blast(qry, subj,window_size,float(json['perc_identity']),gqry,gsubj)
                    
                # Se desejar, você pode chamar a função de plotagem aqui
                # plot_blast_hits(blast_results)
            j += 1
        i += 1


if __name__ == "__main__":
    main()


# In[ ]:




