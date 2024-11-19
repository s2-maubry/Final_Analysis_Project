# Final_Analysis_Project
Extra credit final repository

## Repository Structure
- **Final_Alignments/**: Contains the final alignment file(s) generated during the analysis.
- **Final_Trees/**: Contains the final phylogenetic tree(s) in Newick and PDF formats.
- **Final_Reconciliations/**: Contains the reconciled gene tree(s) and visual representations in SVG and PDF formats.


# Lab 3: Finding homologs with BLAST KEY

Step 1: make a new directory called MARF1
```
mkdir ~/lab03-$MYGIT/METTL2A
```
Step 2: Enter the new directory using the cd command and ensure you are in the right directory
```
cd METTL2A
pwd
```
Step 3: Downlaod METTL2A's protein seuqence in fasta format from NCBI
```
ncbi-acc-download -F fasta -m protein "NP_859076.3"
```
Step 4: Ensure the dowload was successful, you should see NP_859076.3.fa in the METTL2A folder. Next open the file using the less command.
```
ls 
less NP_859076.3.fa
```
Step 5: Perform a BLAST search of METTL2A against the other animal proteomes. The command specifies the output folder and it name.
```
blastp -db ../allprotein.fas -query NP_859076.3.fa -outfmt 0 -max_hsps 1 -out METTL2A.blastp.typical.out
```
Step 6: Look at the output file METTL2A.blastp.typical.out which contains the BLAST results.
```
less METTL2A.blastp.typical.out - look at output file with BLAST results
```
Step 7: Request tabular format for BLAST results and place the new results in a specified output file named METTL2A.blastp.detail.out
```
blastp -db ../allprotein.fas -query NP_859076.3.fa  -outfmt "6 sseqid pident length mismatch gapopen evalue bitscore pident stitle"  -max_hsps 1 -out METTL2A.blastp.detail.out 
```
Step 8: Open the detailed output file containing BLAST results sorted in tabular format
```
less -S MARF1.blastp.detail.out
```
Step 9: Sort the detailed output file to determine how many hits the METTL2A proteome sequence has
```
grep -c H.sapiens MARF1.blastp.detail.out
```
Step 10: We will filter out the hits that have an e-value that is greater then 1e-30, the remaining hits will be placed in a file called METTL2A.blastp.detail.filtered.out. Next, we will count the number of hits remaining in the filtered folder
```
awk '{if ($6< 1e-30)print $1 }' METTL2A.blastp.detail.out > METTL2A.blastp.detail.filtered.out
wc -l MARF1.blastp.detail.filtered.out 
```
> I had 33 homologs with an e-value less than 1e-30

Step 11: Determine how many paralogs were found pers species.
```
grep -o -E "^[A-Z]\.[a-z]+" MARF1.blastp.detail.filtered.out  | sort | uniq -c 
```

# Lab 4: Gene family sequence alignment

Step 1: make a new directory called MARF1
```
mkdir ~/lab04-$MYGIT/METTL2A
```
Step 2: Enter the new directory using the cd command and ensure you are in the right directory
```
cd METTL2A
pwd
```
Step 3: Identify and isolate homologous sequences of interest based on the BLAST results and a specific filtering criteria.
```
seqkit grep --pattern-file ~/lab03-$MYGIT/METTL2A/METTL2A.blastp.detail.filtered.out ~/lab03-$MYGIT/allprotein.fas | seqkit grep -v -p "carpio" > ~/lab04-$MYGIT/METTL2A/METTL2A.homologs.fas       
```
> this loaded 33 patterns from the file

Step 4: Conduct a multiple sequence alignment with the homolog sequences and put this alignment results in the specified output file
```
muscle -align ~/lab04-$MYGIT/METTL2A/METTL2A.homologs.fas -output ~/lab04-$MYGIT/METTL2A/METTL2A.homologs.al.fas
```

Step 5: Convert information into a scrollable format 
```
alv -kli  ~/lab04-$MYGIT/METTL2A/METTL2A.homologs.al.fas | less -RS
```
Step 6: Change seeting to color columns where the most common amino id is in >50% of sequences
```
alv -kli --majority ~/lab04-$MYGIT/METTL2A/METTL2A.homologs.al.fas | less -RS
```
Step 7: Run plot MSA.R using R with the specified alignment file, the results whivh be in a graphical representation of MSA (pdf)
```
Rscript --vanilla ~/lab04-$MYGIT/plotMSA.R  ~/lab04-$MYGIT/METTL2A/METTL2A.homologs.al.fas
```
Step 8: Calculate the width of the alignment 
```
alignbuddy  -al  ~/lab04-$MYGIT/METTL2A/METTL2A.homologs.al.fas      
```
> length of alignment: 793

Step 9: Calculate the length of the alignment after taking out columns with gaps
```
alignbuddy -trm all  ~/lab04-$MYGIT/METTL2A/METTL2A.homologs.al.fas | alignbuddy  -al        
```
> 256 after gaps= 793-256= 537 gaps

Step 10: Calculates the length of the alignment after taking out invariable positions
```
alignbuddy -dinv 'ambig' ~/lab04-$MYGIT/METTL2A/METTL2A.homologs.al.fas | alignbuddy  -al
```
> 741 after invariant positions are removed 793-741= 52

Step 11: T_coffee will calculate the percent indentity excluding gaps
```
t_coffee -other_pg seq_reformat -in ~/lab04-$MYGIT/METTL2A/METTL2A.homologs.al.fas -output sim  
```
> percent indentity using t_coffee was 50.85

Step 12: Use Align buddy to calculate the percent indentity excluding gaps
```
alignbuddy -pi ~/lab04-$MYGIT/METTL2A/METTL2A.homologs.al.fas | awk ' (NR>2)  { for (i=2;i<=NF  ;i++){ sum+=$i;num++} }
     END{ print(100*sum/num) } '
```
> percent indentity using align buddy was 40.5752


# Lab 5: Gene Family Phylogeny using IQ-TREE

Step 1: Make new directory and cd into it
```
mkdir ~/lab05-$MYGIT/METTL2A
cd ~/lab05-$MYGIT/METTL2A
```
Step 2: Remove sequences that contain a duplicate labet tag and make a copy from previous METTL2A folder to this one.
```
sed 's/ /_/g'  ~/lab04-$MYGIT/METTL2A/METTL2A.homologs.al.fas | seqkit grep -v -r -p "dupelabel" >  ~/lab05-$MYGIT/METTL2A/METTL2A.homologsf.al.fas
```
Step 3: Create a maximum likehood tree estimate using IQ-TREE. This comman will first calculate the optimal amino acids substitution and amino acid frequency. Next, it will do a tree search, estimate branch lengths, and estimate ultrafast bootstrap support.
```
iqtree -s ~/lab05-$MYGIT/METTL2A/METTL2A.homologsf.al.fas -bb 1000 -nt 2
```
Step 4: View the .treefile using newick display program
```
nw_display ~/lab05-$MYGIT/METTL2A/METTL2A.homologsf.al.fas.treefile
```
Step 5: Use R to visual the tree in graphical format
```
Rscript --vanilla ~/lab05-$MYGIT/plotUnrooted.R  ~/lab05-$MYGIT/METTL2A/METTL2A.homologsf.al.fas.treefile ~/lab05-$MYGIT/METTL2A/METTL2A.homologsf.al.fas.treefile.pdf 0.4 15
gotree reroot midpoint -i ~/lab05-$MYGIT/METTL2A/METTL2A.homologsf.al.fas.treefile -o ~/lab05-$MYGIT/METTL2A/METTL2A.homologsf.al.mid.treefile
```
Step 6: Use the gotree command to midpoint root the tree to interpret it in a more meaningful way
```
gotree reroot midpoint -i ~/lab05-$MYGIT/METTL2A/METTL2A.homologsf.al.fas.treefile -o ~/lab05-$MYGIT/METTL2A/METTL2A.homologsf.al.mid.treefile
```
Step 7: View the midpoint rooted tree, first on the command line and then in graphical format
```
nw_order -c n ~/lab05-$MYGIT/METTL2A/METTL2A.homologsf.al.mid.treefile  | nw_display -
nw_order -c n ~/lab05-$MYGIT/METTL2A/METTL2A.homologsf.al.mid.treefile | nw_display -w 1000 -b 'opacity:0' -s  >  ~/lab05-$MYGIT/METTL2A/METTL2A.homologsf.al.mid.treefile.svg -
```
Step 8: Convert the midpoint rooted tree in svg format into pdf format 
```
convert  ~/lab05-$MYGIT/METTL2A/METTL2A.homologsf.al.mid.treefile.svg  ~/lab05-$MYGIT/METTL2A/METTL2A.homologsf.al.mid.treefile.pdf
```
Step 9: Switch the view from a phylogram to a cladogram
```
nw_order -c n ~/lab05-$MYGIT/METTL2A/METTL2A.homologsf.al.mid.treefile | nw_topology - | nw_display -s  -w 1000 > ~/lab05-$MYGIT/METTL2A/METTL2A.homologsf.al.midCl.treefile.svg -
convert ~/lab05-$MYGIT/METTL2A/METTL2A.homologsf.al.midCl.treefile.svg ~/lab05-$MYGIT/METTL2A/METTL2A.homologsf.al.midCl.treefile.pdf
```
Step 10: Reroot the tree to include an outgroup
```
nw_reroot ~/lab05-$MYGIT/METTL2A/METTL2A.homologsf.al.fas.treefile H.sapiens_HBG1_hemoglobin_subunit_gamma1 H.sapiens_HBG2_hemoglobin_subunit_gamma2 H.sapiens_HBB_hemoglobin_subunit_beta H.sapiens_HBD_hemoglobin_subunit_delta >~/lab05-$MYGIT/METTL2A/METTL2A.homologsf.outgroupbeta.treefile
```
Step 11: Visualize the new outgroup rooted tree as a pdf
```
nw_order -c n ~/lab05-$MYGIT/METTL2A/METTL2A.homologsf.outgroupbeta.treefile | nw_topology - | nw_display -s  -w 1000 > ~/lab05-$MYGIT/METTL2A/METTL2A.homologsf.outgroupbeta.treefile.svg -

convert ~/lab05-$MYGIT/METTL2A/METTL2A.homologsf.outgroupbeta.treefile.svg ~/lab05-$MYGIT/METTL2A/METTL2A.homologsf.outgroupbeta.treefile.pdf
```

# Lab 6: Reconciling a Gene and Species Tree

Step 1: Make METTL2A directory and cd into it
```
mkdir ~/lab06-$MYGIT/METTL2A
cd METTL2A
```
Step 2: Copy the gene tree from the previous lab into the current directory
```
cp ~/lab05-$MYGIT/METTL2A/METTL2A.homologsf.al.mid.treefile ~/lab06-$MYGIT/METTL2A/METTL2A.homologsf.al.mid.treefile
```

Step 3: Uses java to run NOTUNG and to reconcile the specified gene tree and species tree and to save the results as a png in the specified location, then view the tree
```
java -jar ~/tools/Notung-3.0_24-beta/Notung-3.0_24-beta.jar -s ~/lab05-$MYGIT/species.tre -g ~/lab06-$MYGIT/METTL2A/METTL2A.homologsf.al.mid.treefile --reconcile --speciestag prefix --savepng --events --outputdir ~/lab06-$MYGIT/METTL2A/

less METTL2A.homologsf.al.mid.treefile.rec.events.txt
```

Step 4: Use python to convert the NOTUNG file to RecphyloXML ultimately creating a .xml file
```
python2.7 ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g ~/lab06-$MYGIT/METTL2A/METTL2A.homologsf.al.mid.treefile.rec.ntg --include.species
```

Step 5: Use thirdkind to view the gene tree reconciliation within the species tree, this will display internal nodes, gene duplications, and loss events in an svg file
```
thirdkind -Iie -D 40 -f ~/lab06-$MYGIT/METTL2A/METTL2A.homologsf.al.mid.treefile.rec.ntg.xml -o  ~/lab06-$MYGIT/METTL2A/METTL2A.homologsf.al.mid.treefile.rec.svg
```
Step 6: Convert the svg file of the reconciliation and convert it into a pdf 
```
convert  -density 150 ~/lab06-$MYGIT/METTL2A/METTL2A.homologsf.al.mid.treefile.rec.svg ~/lab06-$MYGIT/METTL2A/METTL2A.homologsf.al.mid.treefile.rec.pdf
```


# Lab 8: Protein Domain Prediction

Step 1: 
```
mkdir ~/lab08-$MYGIT/METTL2A
cd ~/lab08-$MYGIT/METTL2A
```
Step 2: Make a copy of the raw unaligned sequences (with removed stop codons), move the output into a specified folder 
```
sed 's/*//' ~/lab04-$MYGIT/METTL2A/METTL2A.homologs.fas > ~/lab08-$MYGIT/METTL2A/METTL2A.homologs.fas
```
Step 3: Run the output file through rpsblast, we will do this by first copying the unligned sequence then use rpsblast to compare the query protein to a database of pre computed profiles; this will help identify protein domains. The output file will give information about the length, start, and end of alignmnet as well as give an e-value to determine the significance of the matches.
```
rpsblast -query ~/lab08-$MYGIT/METTL2A/METTL2A.homologs.fas -db ~/data/Pfam/Pfam -out ~/lab08-$MYGIT/METTL2A/METTL2A.rps-blast.out  -outfmt "6 qseqid qlen qstart qend evalue stitle" -evalue .0000000001
```

Step 4: Copy final gene tree from a previous lab into the current folder.
```
cp ~/lab05-$MYGIT/METTL2A/METTL2A.homologsf.outgroupbeta.treefile ~/lab08-$MYGIT/METTL2A
```
Step 5: Use R to determine predicted domains for the query protein, in order to do this 3 files must be specified: midpoint rooted tree file, predicted pfram domains from rpsblast, and the specification for the name of the output pdf we want to write.
```
Rscript  --vanilla ~/lab08-$MYGIT/plotTreeAndDomains.r ~/lab08-$MYGIT/METTL2A/METTL2A.homologsf.outgroupbeta.treefile ~/lab08-$MYGIT/METTL2A/METTL2A.rps-blast.out ~/lab08-$MYGIT/METTL2A/METTL2A.tree.rps.pdf
```

Step 6: Convenient way to look at the annotations located in the rps-blast.out file in a spreadsheet program.
```
mlr --inidx --ifs "\t" --opprint  cat ~/lab08-$MYGIT/METTL2A/METTL2A.rps-blast.out | tail -n +2 | less -S
```

Step 7: Count how many proteins have more than one annotation, determine which pfam domain annotations is found most often, and which proteins has the longest annotated protein domain
Which Pfam domain annotation is most commonly found
```
cut -f 1 ~/lab08-$MYGIT/METTL2A/METTL2A.rps-blast.out | sort | uniq -c

cut -f 6 ~/lab08-$MYGIT/METTL2A/METTL2A.rps-blast.out | sort | uniq -c

awk '{a=$4-$3;print $1,'\t',a;}' ~/lab08-$MYGIT/METTL2A/METTL2A.rps-blast.out |  sort  -k2nr
```
#### Most amount of pfam domains
> 8 D.rerio_mettl2a_tRNA_N3methylcytidine_methyltransferase_METTL2, this gene in D.rerio has 8 pfam domains.

#### Pfam domains found in the most genes
> 18 pfam08241, Methyltransf_11, Methyltransferase domain.  Members of this family are SAM dependent methyltransferases.
> 34 pfam08242, Methyltransf_12, Methyltransferase domain.  Members of this family are SAM dependent methyltransferases.
> 13 pfam13489, Methyltransf_23, Methyltransferase domain.  This family appears to be a methyltransferase domain.
> 34 pfam13649, Methyltransf_25, Methyltransferase domain.  This family appears to be a methyltransferase domain.

#### Pfam domains found multiple times in a single gene
> The Methyltransferase family has various domains specifically, Methyltransf_12 and Methyltransf_25 which are each found 34 times, making them common in my analysis. Additionaly, the gene D.rerio_mettl2a_tRNA_N3methylcytidine_methyltransferase_METTL2 contains multiple Pfam domains which are pfam08241, pfam08242, pfam13489, and pfam13649. These each appear twice.

#### Gene with the longest pfam domain annotation
> C.mydas_METTL6_tRNA_N3methylcytidine_methyltransferase_METTL6_iso  177. This gens has the longest pfam domain, spanning 177 residues. This domain is pfam13489 (found this by going to the rps-blast.out file and finding C.mydas, then since there were 3 looking which one had the longest length)

#### Domain annotation with the best e-value
> C.carcharias_mettl6_tRNA_N3methylcytidine_methyltransferase_METTL has an e-value of 1.55e-21. This is the best E-value because it is the smallest, indicating it has the lowest chance of being an alignment due to random chance.

Step 8: Continuing analysis to determine which protein has the shortest annotated protein domain, and determine hich protein has a domain annotation with the best e-value
```
cut -f 1,5 -d $'\t' ~/lab08-$MYGIT/METTL2A/METTL2A.rps-blast.out | sort -t $'\t' -k2,2g

cut -f 1,6 ~/lab08-$MYGIT/METTL2A/METTL2A.rps-blast.out | tr ',' '\n' | sort | uniq -c | awk '$1 > 1'

cut -f 1,6 ~/lab08-$MYGIT/METTL2A/METTL2A.rps-blast.out | awk '{print length($2), $0}' | sort -nr | head -n 1

cut -f 1,6 ~/lab08-$MYGIT/METTL2A/METTL2A.rps-blast.out | sort | uniq -c
```








