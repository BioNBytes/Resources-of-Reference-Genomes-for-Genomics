🧬 Reference Genomes and Genome Annotations: Guidelines for Human and Mouse
A curated guide for selecting, understanding, and using reference genomes and gene annotations in genomic studies.

📌 Overview
Reference genomes and genome annotations are foundational to all genomic analyses. From RNA-seq to variant calling, the choice of genome version, annotation format, and source database directly influences the quality, reproducibility, and interpretability of results.

This guide aims to:

Introduce concepts of genome assembly and annotation.

Compare key platforms (NCBI, UCSC, Ensembl, GENCODE).

Provide access examples for retrieving mouse gene data.

Discuss practical considerations when choosing annotations for analysis.

🧱 Basics of Genome Assembly
From Reads to Chromosomes
Contigs: Continuous sequences assembled from reads.

Scaffolds: Ordered and oriented contigs.

Chromosomes: Final level of assembly (may include placed/unlocalized/unplaced scaffolds).

Scaffold Classification
Type	Description
Placed	Positioned within a chromosome.
Unlocalized	Assigned to a chromosome but location uncertain.
Unplaced	Chromosome of origin unknown.

📖 Learn more at Ensembl

🧬 Reference Genome Versions
Human (Homo sapiens)
GRCh37: Widely used legacy version; includes alternate loci and patches.

GRCh38: Current standard reference; includes 178 alternate loci (e.g., MHC).

T2T-CHM13: Telomere-to-Telomere complete genome (lacks Y chromosome).

Pangenome: Ongoing project capturing human diversity.

Mouse (Mus musculus)
GRCm39 (mm39): Current reference genome used in most modern pipelines.

🧾 Genome Annotation Formats
File Types
Format	Notes
GTF (Gene Transfer Format)	Tab-delimited, simpler, often used with RNA-seq pipelines.
GFF3 (General Feature Format v3)	More structured, supports multiple levels of annotation.

Key Challenges
Ambiguity in overlapping genes.

Format inconsistencies (especially UCSC GTFs).

Differences in transcript biotypes and gene models across databases.

🔍 Major Platforms: Comparison & Use Cases
Feature	NCBI (RefSeq)	UCSC	Ensembl	GENCODE
Species	Broad	Broad	Broad	Human, Mouse only
Annotation Type	Curated + automated	Known genes (protein-focused)	Evidence-based + comparative	Comprehensive, harmonized with Ensembl
Strengths	Official RefSeq IDs, good NCBI integration	Fast browser, multiple data types	Rich comparative genomics	Gold standard for transcript annotations
Limitations	Fewer transcript isoforms	GTF inconsistencies	Complex structure for newcomers	Only human and mouse

🧠 “No database is perfect. Some gene annotations may be inaccurate or even wrong.”

⚙️ Genome Annotation & Quantification Tips
For RNA-seq:

Use simpler annotations for robust quantification (e.g., remove non-expressed or redundant genes).

Prefer GENCODE for human/mouse expression analysis.

For novel transcript discovery:

Use Ensembl or GENCODE full annotation.

For protein-coding focus:

UCSC Known Genes or NCBI RefSeq may be appropriate.

🧪 Impact of annotation complexity on differential expression analysis (PMID: 34625021)

🧰 Accessing Genome Data Programmatically
🔗 NCBI Entrez Direct: Mouse Gene Access
bash
Copy
Edit
# Get mouse gene info (e.g., TP53)
esearch -db gene -query "Tp53[Gene Name] AND Mus musculus[Organism]" | \
elink -target nuccore | \
efetch -format fasta
🧬 Ensembl REST API: Mouse Gene Info
python
Copy
Edit
import requests

gene = "Tp53"
server = "https://rest.ensembl.org"
ext = f"/lookup/symbol/mus_musculus/{gene}?expand=1"

r = requests.get(server+ext, headers={ "Content-Type" : "application/json" })
if not r.ok:
    r.raise_for_status()

decoded = r.json()
print(f"Gene ID: {decoded['id']}")
print(f"Transcript IDs: {[t['id'] for t in decoded.get('Transcript', [])]}")
📚 Additional Resources
NCBI Genome Database

Ensembl Genome Browser

UCSC Genome Browser

GENCODE Project

Genome Assembly Overview

INSCD Collaboration

RNA-seq Quantification and Annotations (Review)

🧩 Final Thoughts
Choosing the right reference genome and annotation is not trivial—it depends on your analysis goals (e.g., expression quantification vs. variant discovery), species of interest, and available tools.

“Annotation complexity can limit reproducibility, but oversimplification may miss important biology.”
