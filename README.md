## Resources on Reference Genomes and Annotations

This guide contains key concepts and summary of resources on reference genomes and annotations for human and mouse.

### 1. Reference Genome Assembly: Concepts and Hierarchies

**Genome assembly** is the process of reconstructing the complete DNA sequence of an organism's genome from short or long sequencing reads. The hierarchy of assembly follows a structured progression:

| Level                    | Description                                                                                                                                         |
| ------------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Reads**                | Raw sequence data from sequencing platforms (short e.g. Illumina or long reads e.g. PacBio, Oxford Nanopore).                                       |
| **Contigs**              | Continuous sequences formed by overlapping and merging reads. No gaps.                                                                              |
| **Scaffolds**            | Larger, ordered and oriented contigs, separated by estimated gaps (Ns).                                                                             |
| **Chromosomes / Genome** | Full assembly where scaffolds are assigned to chromosomes, created by anchoring scaffolds using additional data (Hi-C, optical maps, linkage maps). |

**Types of Scaffolds:**

|Scaffold Type|Description|
|---|---|
|**Placed scaffolds**|Assigned to a precise chromosomal location.|
|**Unlocalized scaffolds**|Assigned to a specific chromosome, but position/orientation is uncertain.|
|**Unplaced scaffolds**|Cannot be confidently assigned to any chromosome. Often contain repeats or unresolved sequences.|

**Common Genome Assembly Related Terms and their Definitions:**
- **Reference Genome**  
    A curated, representative sequence of a species’ genome, typically derived from one or a few individuals.
- **Top-level sequence**  
    The highest-level components of a genome assembly, including chromosomes, scaffolds, and contigs, used to reconstruct the genome.
- **Primary assembly**  
    The main genome assembly version, excluding alternate haplotypes or loci; includes placed chromosomes and unlocalized/unplaced scaffolds.
- **Alternate loci**  
    Alternate haplotypes or divergent regions (e.g., MHC), representing genetic variation within the species.
- **Patch sequences**  
    Sequence corrections or updates that do not change the chromosome coordinate system.
### 2. Human Reference Genomes

The **human reference genome** is a digital representation of human DNA used for gene annotation, variant calling, and sequence alignment. The most commonly used assemblies, **GRCh37 (hg19)** and **GRCh38 (hg38)**, are **linear assemblies** derived from a small number of individuals, primarily of European ancestry. GRCh38 added modeled centromeres, alternate loci for highly polymorphic regions (like the MHC), and improved gap representation. However, both assemblies have limitations: they lack full coverage of repetitive regions (e.g., centromeres, telomeres), suffer from **reference bias** (favoring alleles present in the reference), and underrepresent global genetic diversity.

Key features include the use of **top-level sequences**—which comprise chromosomes, scaffolds, and unplaced/unlocalized contigs—and integration with annotation frameworks like **GENCODE**, **Ensembl**, and **RefSeq**. **Alternate loci** provide optional representations of highly variable regions. To overcome the shortcomings of linear references, the **T2T-CHM13** assembly (2022) delivered the first **gapless, telomere-to-telomere** genome, though it is haploid and lacks heterozygosity. The **Human Pangenome Reference Consortium (HPRC)** is now building a **graph-based reference**, which incorporates hundreds of diverse haplotypes and allows multiple paths through complex genomic regions, thereby reducing mapping bias and better reflecting population diversity.

Despite these advancements, GRCh38 remains the dominant reference in most genomic pipelines, underscoring the need for broader adoption of complete and inclusive references like T2T and pangenomes.

| Assembly Version     | Release      | Features                                                                                                     |
| -------------------- | ------------ | ------------------------------------------------------------------------------------------------------------ |
| **GRCh37 (hg19)**    | 2009         | Widely used, lacks complete centromeres and telomeres. Still used in legacy pipelines.                       |
| **GRCh38 (hg38)**    | 2013         | Improved centromere modeling, added alternate loci (e.g., for MHC), better gap representation.               |
| **T2T-CHM13**        | 2022         | First telomere-to-telomere (T2T) assembly. No gaps, complete centromeres, acrocentric arms, and rDNA arrays. |
| **pangenome (HPRC)** | 2023–present | Graph-based genome capturing population diversity across many haplotypes and ancestries.                     |


### 3. Mouse Reference Genomes
The **mouse reference genome** has undergone several important updates over the years. It began with the initial draft assembly, **MGSCv37**, released in 2004, followed by improvements in **NCBIm37** (2007) and **GRCm38** (2011). The current standard, **GRCm39**, was released in 2020. These assemblies are mainly based on the **C57BL/6J inbred strain**, a widely used laboratory mouse model, and provide high-quality, **chromosome-level sequences**.

The reference genome includes well-annotated **genes** and **regulatory elements**, supported by resources such as **GENCODE**, **Ensembl**, and **RefSeq**. It consists of **top-level sequences** like **chromosomes**, **scaffolds**, and **unplaced/unlocalized contigs**. Some **alternate loci** are included to represent strain-specific variants, though these are limited compared to the diversity seen across mouse strains.

Despite its utility, the mouse reference genome has **limitations**. Because it represents a single inbred strain, it does not capture the broad **genetic diversity** found in other laboratory and wild-derived mice. This results in **reference bias** during read mapping and variant calling when analyzing other strains. Additionally, complex regions such as **repetitive sequences**, **centromeres**, and **structural variants** remain incompletely resolved.

To overcome these challenges, ongoing efforts like the **Mouse Pangenome Project** aim to develop **graph-based**, **multi-strain reference genomes**. These will incorporate sequences from multiple mouse strains to better represent genetic variation, enabling more accurate genomic analyses across diverse mouse populations.

| Assembly Version              | Release Year | Key Features                                                                                               |
| ----------------------------- | ------------ | ---------------------------------------------------------------------------------------------------------- |
| **MGSCv37 (mm5)**             | 2004         | Early draft assembly; fragmented, incomplete chromosomes.                                                  |
| **NCBIm37 (mm9)**             | 2007         | Improved assembly; used extensively in early mouse studies.                                                |
| **GRCm38 (mm10)**             | 2011         | Major update; improved contiguity, added alternate loci; widely used.                                      |
| **GRCm39 (mm39)**             | 2020         | Current standard; chromosome-level assembly based on C57BL/6J strain; improved annotation and gap closure. |
| **Mouse Pangenome (ongoing)** | Ongoing      | Graph-based, multi-strain references under development; aims to capture genetic diversity across strains.  |
|                               |              |                                                                                                            |

### 4. Annotation Resources
The choice of annotation source and format directly impacts mapping, quantification, and interpretation in genomics. Here is a detailed, accurate comparison of the most widely used resources:

| **Resource**            | **Source / Consortium**                          | **Annotation Set Name**                 | **Formats**              | **Transcript/ID Format**                            | **Features**                                                                                                                                                                 | **Limitations**                                                 |
| ----------------------- | ------------------------------------------------ | --------------------------------------- | ------------------------ | --------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------- |
| **NCBI RefSeq**         | NCBI                                             | RefSeq                                  | GFF, GTF, GenBank, FASTA | NM_ (mRNA), NR_ (ncRNA), XP_, NP_ (protein)         | Manually curated and computationally predicted gene models; broad organism coverage; non-redundant, well-annotated; provides a "representative" transcript per gene.         | Limited species breadth compared to Ensembl; slower updates     |
| **Ensembl**             | EMBL-EBI                                         | Ensembl Genes                           | GTF, GFF3, FASTA         | ENSG (gene), ENST (transcript), ENSP (protein)      | Broad species coverage; integrates GENCODE, RefSeq, UniProt, RNA-seq; frequent updates; supports genome browser and API access.                                              | Automated pipelines may miss rare isoforms; update lag          |
| **GENCODE**             | GENCODE Consortium (EMBL-EBI & Sanger Institute) | GENCODE                                 | GTF, GFF3, FASTA         | ENSG/ENST (same as Ensembl)                         | Most comprehensive, manually curated annotation for human and mouse; merges Ensembl and HAVANA; includes protein-coding, lncRNA, pseudogenes, sRNAs; default in GTEx/ENCODE. | Focus on human/mouse only; complex format for beginners         |
| **UCSC Genome Browser** | UCSC                                             | UCSC Genes, refGene, ensGene, knownGene | BED, GTF, BigBed, BAM    | uc#####.# (UCSC), NM_/NR_ (refGene), ENST (ensGene) | Integrates multiple annotation datasets (RefSeq, GENCODE, Ensembl); powerful visualization tools; customizable with user tracks.                                             | Quality varies by source; large data volume can be overwhelming |
| **MANE**                | NCBI + Ensembl                                   | MANE Select / MANE Plus Clinical        | GTF, FASTA               | NM_/ENST (matched)                                  | One matched transcript per protein-coding gene, agreed by RefSeq and Ensembl; ideal for clinical/diagnostic use; high concordance.                                           | Limited to protein-coding genes                                 |
| **UniProt**             | UniProt Consortium                               | UniProt                                 | FASTA, XML               | Protein accessions                                  | Protein sequence and functional annotation database; links genomic data to protein function, structure, and interactions.                                                    | Protein-centric; less focused on transcript structures          |

**Additional Notes and Details:**
- **GENCODE** (EMBL-EBI & Sanger Institute): High-quality, comprehensive annotation for human and mouse combining manual curation and automated pipelines; regularly updated and default for projects like GTEx and ENCODE. Mainly focused on human/mouse; complex formats can challenge beginners.
- **Ensembl** (EMBL-EBI): Broad species coverage with automated, evidence-based annotation integrating GENCODE, RefSeq, UniProt, and RNA-seq; supports genome browsing and API. Automated pipelines may miss rare isoforms; transcript IDs can change.
- **NCBI RefSeq**: Curated, non-redundant gene and transcript models for 165,000+ species; combines manual and computational curation; provides representative transcripts. Smaller species range and slower updates than Ensembl.
- **UCSC Genome Browser** (UC Santa Cruz): Integrates annotation tracks from RefSeq, GENCODE, Ensembl, and comparative genomics with customizable visualization and user track support. Quality varies; transcript IDs differ by track (e.g., uc#####.#, NM_/NR_, ENST).
- **MANE** (NCBI + Ensembl): Collaborative set providing one matched transcript per protein-coding gene for clinical use, ensuring high concordance; limited to protein-coding genes.
- **UniProt**: Protein sequence and functional annotation linking genomic data to protein function and interactions; protein-centric and less detailed on transcript isoforms/non-coding RNAs.

**Annotation Formats:**

| Format   | Description & Use Case                                                                                                                           |
| -------- | ------------------------------------------------------------------------------------------------------------------------------------------------ |
| **GTF**  | Widely used, simpler format; supports gene and transcript features; ideal for RNA-seq quantification and legacy tools.                           |
| **GFF3** | More structured and flexible; supports hierarchical relationships (gene→transcript→exon); better for complex feature annotation and submissions. |
| **BED**  | Simple format for genomic regions; often used for visualization and custom tracks.                                                               |

**Impact on Analysis:** 
- Annotation sets vary in **transcript models**, **exon-intron structures**, and **curation methods**, which affects RNA-seq mapping, quantification, and downstream interpretation.
- Only a minority of genes show identical quantification across resources due to these differences.
- For **human and mouse** transcriptomics, **GENCODE** is generally recommended; **NCBI RefSeq** is favored for clinical or highly curated applications; **Ensembl** suits broad discovery; **UCSC Table Browser** is useful for custom, region-focused analyses.
- Always document the **annotation source**, **version**, and **genome build** used in analyses.
- Format inconsistencies (especially with some UCSC GTFs) may cause downstream pipeline issues.


### 5. Programmatic Access: Mouse Gene Example

Below are code examples for retrieving canonical transcript IDs for mouse or human genes from both NCBI and Ensembl.

#### A. Retrieving Canonical NCBI RefSeq Transcript for a Mouse Gene
- Uses NCBI E-utilities (`esearch`, `elink`, `esummary`) to find curated NM_ transcripts by gene symbol and organism.
- Includes retry with exponential backoff for HTTP 429 (rate limit) errors.
- Works for both **mouse** and **human** genes.**

```python
import time
import requests
from xml.etree import ElementTree as ET

def make_request_with_retry(url, params, max_retries=5):
    """
    Make an HTTP GET request with retries on HTTP 429 (Too Many Requests).
    Uses exponential backoff (wait times: 1, 2, 4, 8, 16 seconds).
    """
    for attempt in range(max_retries):
        response = requests.get(url, params=params)
        if response.status_code == 429:
            wait = 2 ** attempt  # exponential backoff: 1,2,4,8,16 seconds
            print(f"429 received, sleeping for {wait} seconds...")
            time.sleep(wait)
        else:
            return response  # Successful or other error code
    # If max retries reached without success, raise an exception
    raise Exception("Max retries reached with 429 responses")

def get_ncbi_canonical_transcript(gene_symbol, organism="Mus musculus"):
    """
    Retrieve the canonical NCBI RefSeq transcript (NM_ accession) for a given gene symbol and organism.
    Steps:
      1. Search for gene ID using esearch.
      2. Find linked RefSeq RNA transcripts via elink.
      3. Summarize transcript info using esummary.
      4. Return the first curated NM_ transcript found.
    """
    # Step 1: Use esearch to find the gene ID for the gene symbol and organism
    search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    search_params = {
        "db": "gene",
        "term": f"{gene_symbol}[Gene Name] AND {organism}[Organism]",
        "retmode": "json"
    }
    response = make_request_with_retry(search_url, search_params)
    
    if not response.ok:
        return f"NCBI esearch request failed: {response.status_code} - {response.reason}"

    try:
        json_data = response.json()
    except ValueError:
        return f"Failed to parse JSON from NCBI esearch response:\n{response.text}"

    esearch_result = json_data.get("esearchresult")
    if not esearch_result:
        return f"No 'esearchresult' found in response:\n{json_data}"

    gene_id_list = esearch_result.get("idlist", [])
    if not gene_id_list:
        return f"No gene found for {gene_symbol} in {organism}"

    gene_id = gene_id_list[0]

    # Step 2: Use elink to find RefSeq RNA transcripts linked to this gene ID
    link_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
    link_params = {
        "dbfrom": "gene",
        "db": "nucleotide",
        "id": gene_id,
        "retmode": "json",
        "linkname": "gene_nuccore_refseqrna"
    }
    response = make_request_with_retry(link_url, link_params)

    try:
        linked_ids = response.json()["linksets"][0]["linksetdbs"][0]["links"]
    except (KeyError, IndexError):
        return f"No RefSeq transcripts found for {gene_symbol} in {organism}"

    if not linked_ids:
        return f"No linked RefSeq transcripts found for {gene_symbol} in {organism}"

    # Step 3: Use esummary to get transcript details for linked transcripts
    ids_str = ",".join(linked_ids[:50])  # Limit to first 50 transcripts
    summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    summary_params = {
        "db": "nucleotide",
        "id": ids_str,
        "retmode": "xml"
    }
    response = make_request_with_retry(summary_url, summary_params)

    root = ET.fromstring(response.content)
    transcripts = []
    for docsum in root.findall(".//DocSum"):
        accession = None
        for item in docsum.findall("Item"):
            if item.attrib.get("Name") == "Caption" and item.text.startswith("NM_"):
                accession = item.text
        if accession:
            transcripts.append(accession)

    # Step 4: Return the first curated NM_ transcript found or an error message
    if transcripts:
        return transcripts[0]
    else:
        return f"No NM_ (curated) transcripts found for {gene_symbol} in {organism}"


# Example usage for mouse
mouse_gene = "Kras"
mouse_organism = "Mus musculus"
mouse_canonical = get_ncbi_canonical_transcript(mouse_gene, mouse_organism)
print(f"Canonical NCBI transcript for {mouse_gene} ({mouse_organism}): {mouse_canonical}")


# Example usage for human
human_gene = "KRAS"
human_organism = "Homo sapiens"
human_canonical = get_ncbi_canonical_transcript(human_gene, human_organism)
print(f"Canonical NCBI transcript for {human_gene} ({human_organism}): {human_canonical}")
```

#### **B. Retrieving Canonical Ensembl Transcript for a Mouse Gene**
- Uses Ensembl REST API to get canonical transcript by gene symbol and species (`mus_musculus` or `homo_sapiens`).
- Returns Ensembl canonical transcript ID.
- Supports both mouse and human.

```python
import requests

def get_ensembl_canonical_transcript(gene_symbol, species="mus_musculus"):
    """
    Retrieve the canonical Ensembl transcript ID for a given gene symbol and species.

    Parameters:
        gene_symbol (str): Gene symbol (e.g., "Kras", "BRCA1").
        species (str): Species name in Ensembl format (default: "mus_musculus" for mouse).
                       For human use "homo_sapiens".

    Returns:
        str: Canonical transcript ID or error message if not found.
    """
    server = "https://rest.ensembl.org"
    ext = f"/lookup/symbol/{species}/{gene_symbol}?expand=1"

    headers = {"Content-Type": "application/json"}
    response = requests.get(server + ext, headers=headers)

    if not response.ok:
        return f"Error fetching data: {response.status_code} - {response.reason}"

    data = response.json()

    canonical_transcript_id = data.get("canonical_transcript")
    if canonical_transcript_id:
        return canonical_transcript_id
    else:
        return f"No canonical transcript found for {gene_symbol} in {species}"


# Example usage for mouse
mouse_gene = "Kras"
mouse_organism = "Mus musculus"
mouse_canonical = get_ensembl_canonical_transcript(mouse_gene, mouse_organism)
print(f"Canonical NCBI transcript for {mouse_gene} ({mouse_organism}): {mouse_canonical}")

# Example usage for human
human_gene = "KRAS"
human_organism = "Homo sapiens"
human_canonical = get_ensembl_canonical_transcript(human_gene, human_organism)
print(f"Canonical NCBI transcript for {human_gene} ({human_organism}): {human_canonical}")

```

**Additional Tips for Programmatic Access:**
- **API Rate Limits:** Both NCBI and Ensembl APIs have request limits. For large-scale queries, there maybe delays.
- **Canonical Transcript Definitions:** The definition of "canonical" may differ between databases
- **Batch Processing:** Both scripts can be adapted for batch processing by iterating over a list of gene symbols.

### 6. Additional Key Points and Resources
- **Species-specific databases** (e.g., Wormbase, Flybase) and AWS iGenomes exist, but may be outdated.
- **Chromosome naming conventions** differ by platform (e.g., NCBI/Ensembl: 1,2,...,X,Y,MT; UCSC: chr1, chr2, ..., chrM).
- **Annotation files should use stable gene IDs** (not gene symbols) for consistency across versions.
- **Quality control**: Use annotation files with a `gene_biotype` field for best results in automated pipelines.
- **Download locations**: NCBI, Ensembl, GENCODE, and UCSC all provide FTP and browser-based access to assemblies and annotations.

### References
1. Haibo Liu, LinkedIn. [Notes on reference genome resources](https://www.linkedin.com/feed/update/urn:li:activity:7309376644169752577/?originTrackingId=1OaGKrcAQMub5NqA11hTZw%3D%3D)
2. UCSC News: Human pangenome. https://news.ucsc.edu/2023/05/pangenome-draft/
3. NCBI Insights: GRCm39. https://ncbiinsights.ncbi.nlm.nih.gov/2020/10/16/refseq-grcm39/
4. Comparison of human (and other) genome browsers. PMC3525149. https://pmc.ncbi.nlm.nih.gov/articles/PMC3525149/
5. excluderanges: exclusion sets for T2T-CHM13, GRCm39. https://www.biorxiv.org/content/10.1101/2022.11.21.517407v1.full.pdf
6. GENCODE 2025: reference gene annotation for human and mouse. Nucleic Acids Res 2025;53:D966-D975. https://pubmed.ncbi.nlm.nih.gov/39565199/
7. Morales J, et al. A joint NCBI and EMBL-EBI transcript set for clinical genomics and research. Nature 2022;604:310–315.
8. GENCODE Publications. https://www.gencodegenes.org/pages/publications.html
9. NCBI: Annotating Genomes with GFF3 or GTF files. https://www.ncbi.nlm.nih.gov/genbank/genomes_gff
10. A comprehensive evaluation of Ensembl, RefSeq, and UCSC annotations in the context of RNA-seq read mapping and gene quantification. BMC Genomics. 2015;16:97.
11. Feature annotation: RefSeq vs Ensembl vs Gencode, what's the difference? Bioinformatics Stack Exchange, 2017.
12. The UCSC Table Browser data retrieval tool. Nucleic Acids Res. 2004;32(Database issue):D493–D496.
13. UCSC Genome Browser (and refGene) vs NCBI Gene (RefSeq). UCSC Genome Browser Public Support, 2013.
