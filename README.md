## Reference Genomes and Genome Annotations

This version incorporates and verifies information from the provided markdown, the attached PDF, and authoritative sources. It corrects, clarifies, and expands on key points, especially regarding assembly versions, annotation formats, database strengths/limitations, and best practices.

### **1. Reference Genome Assembly: Concepts and Hierarchies**

- **Genome assembly** proceeds from sequencing reads → contigs → scaffolds → chromosomes, using technologies like BioNano optical mapping and Hi-C for scaffolding[^1].
- **Scaffold types**:
    - **Placed**: Located at a specific position within a chromosome.
    - **Unlocalized**: Assigned to a chromosome but with uncertain position/orientation.
    - **Unplaced**: Chromosome of origin unknown[^1].
- **Toplevel sequences** in Ensembl are all primary chromosomes plus any unlocalized/unplaced scaffolds[^1].


### **2. Human Reference Genomes**

| Version | Release | Key Features \& Limitations |
| :-- | :-- | :-- |
| **GRCh37** | 2009 | Legacy, includes alternate loci and patches, still used for compatibility. |
| **GRCh38** | 2013 | Current standard, 178 alternate loci (e.g., MHC), improved coverage, reduced errors[^1][^2][^3]. |
| **T2T-CHM13** | 2022 | First telomere-to-telomere assembly (except Y), filled last 8% of genome, added ~200Mbp, 1956 new genes[^4]. |
| **Pangenome** | 2023+ | Graph-based, incorporates 47+ individuals, adds 119Mbp, enhances variant detection and diversity[^2][^3]. |

**Limitations**:

- T2T-CHM13 lacks the Y chromosome and, while more complete, represents a single haplotype, not population diversity[^4][^2].
- The pangenome addresses reference bias by representing a broader spectrum of human variation[^2][^3].


### **3. Mouse Reference Genomes**

| Version | Release | UCSC Name | Notes |
| :-- | :-- | :-- | :-- |
| **GRCm39** | 2020 | mm39 | Current, major coordinate update, improved N50, closed gaps, new gene models, based on C57BL/6J[^5][^6][^7]. |
| **GRCm38** | 2011 | mm10 | Previous standard, still widely used in legacy datasets[^8][^5][^6][^7]. |

- GRCm39 resolved 400+ issues, nearly doubled scaffold N50, and added 1.9 Mb of sequence[^5].
- NCBI RefSeq Annotation Release 109 for GRCm39 includes >45,000 curated transcripts and new gene predictions[^5].


### **4. Major Platforms: Roles and Comparison**

| Platform | Species Support | Annotation Basis | Strengths | Limitations |
| :-- | :-- | :-- | :-- | :-- |
| **NCBI (RefSeq)** | Broad | Curated + automated | Integration, stability, curation | Fewer isoforms, slower updates |
| **UCSC** | Broad | Protein-centric, known genes | Visualization, custom tracks | GTF inconsistencies, protein focus |
| **Ensembl** | Broad | Evidence + comparative | Broad transcript set, APIs | Complex for newcomers |
| **GENCODE** | Human, Mouse | Ensembl + manual curation | Gold-standard for transcripts | Only two species |

**Note:** No database is perfect; gene models may differ and some annotations are inaccurate or incomplete[^1][^9][^10].

### **5. Genome Annotation Formats: GTF vs GFF3**

| Format | Features | Best Use Cases |
| :-- | :-- | :-- |
| **GTF** | Simpler, widely used for RNA-seq, two-level hierarchy (gene→transcript), flexible[^1][^11]. | Quantification, legacy tools |
| **GFF3** | Structured, supports multi-level relationships (gene→transcript→exon), richer metadata[^1][^9][^11]. | Submissions, complex features |

- GFF3 is more standardized and supports richer annotation, but both are accepted for many pipelines[^9][^11].
- Format inconsistencies (especially in UCSC GTFs) can cause downstream issues[^1][^9][^10].


### **6. Annotation Choice: Impact and Best Practices**

- **Annotation source and complexity** directly affect RNA-seq mapping, quantification, and differential expression[^1][^10].
- **GENCODE** is recommended for human/mouse expression studies due to its quality and completeness[^10].
- For robust quantification, use simpler annotations and exclude non-mRNA features (miRNA, rRNA, tRNA) unless needed[^1][^10].
- Remove unexpressed gene models and use stranded RNA-seq for better assignment[^1].
- For exploratory work or novel transcript discovery, more complex annotations may be warranted, but this can reduce statistical power in differential analysis[^1][^10].


### **7. Special Topics**

**Alternate Loci and Patches (Human):**

- GRCh38 includes 178 alternate loci for regions like MHC, allowing representation of population variation[^1].
- Patches (novel/fix) provide updates without changing coordinate systems; incorporated into next major releases[^1].

**Decoy Genomes:**

- Used to reduce misalignments, especially for repetitive regions or unassembled sequence[^1].

**Graph-based Pangenomes:**

- Represent multiple haplotypes/alleles, improving variant detection and reducing reference bias[^2][^3].


### **8. Programmatic Access: Mouse Gene Example**

Below are code examples for retrieving various types of information (e.g., canonical transcript IDs) for mouse genes (can be adapted to be used for human genome as well) from both NCBI and Ensembl.


#### **A. Retrieving Canonical NCBI RefSeq Transcript for a Mouse Gene**

This Python script uses the NCBI E-utilities API to:

- Search for the mouse gene symbol and retrieve its GeneID.
- Find all RefSeq RNA transcripts linked to that gene.
- Filter for curated NM_ RefSeq mRNA accessions.
- Return the first NM_ accession as the canonical transcript (commonly used convention).

```python
import requests
from xml.etree import ElementTree as ET

def get_ncbi_canonical_transcript(gene_symbol, organism="Mus musculus"):
    # Step 1: Use esearch to get the Gene ID for the given gene symbol and organism
    search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    search_params = {
        "db": "gene",
        "term": f"{gene_symbol}[Gene Name] AND {organism}[Organism]",
        "retmode": "json"
    }
    response = requests.get(search_url, params=search_params)
    gene_id_list = response.json()["esearchresult"]["idlist"]

    if not gene_id_list:
        return f"No gene found for {gene_symbol}"

    gene_id = gene_id_list[0]

    # Step 2: Use elink to find RefSeq RNA nucleotide records linked to the Gene ID
    link_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
    link_params = {
        "dbfrom": "gene",
        "db": "nucleotide",
        "id": gene_id,
        "retmode": "json",
        "linkname": "gene_nuccore_refseqrna"
    }
    response = requests.get(link_url, params=link_params)
    try:
        linked_ids = response.json()["linksets"][0]["linksetdbs"][0]["links"]
    except (KeyError, IndexError):
        return f"No RefSeq transcripts found for {gene_symbol}"

    # Step 3: Use esummary to get detailed transcript info (limit to first 50 IDs)
    ids_str = ",".join(linked_ids[:50])
    summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    summary_params = {
        "db": "nucleotide",
        "id": ids_str,
        "retmode": "xml"
    }
    response = requests.get(summary_url, params=summary_params)
    root = ET.fromstring(response.content)

    # Step 4: Filter for curated NM_ (RefSeq mRNA) accessions
    transcripts = []
    for docsum in root.findall(".//DocSum"):
        accession = None
        for item in docsum.findall("Item"):
            if item.attrib.get("Name") == "Caption":
                if item.text.startswith("NM_"):
                    accession = item.text
        if accession:
            transcripts.append(accession)

    if transcripts:
        # Return the first NM_ accession as canonical (commonly used convention)
        return transcripts[0]
    else:
        return f"No NM_ (curated) transcripts found for {gene_symbol}"

# Example usage
gene_symbol = "Kras"
canonical_transcript_id = get_ncbi_canonical_transcript(gene_symbol)
print(f"Canonical NCBI transcript for {gene_symbol}: {canonical_transcript_id}")
```

**Notes:**

- NM_ accessions represent curated mRNA transcripts in RefSeq and are widely used as canonical references.
- For some genes, multiple NM_ accessions may exist; the first is often used by convention, but users should check NCBI for gene-specific canonical definitions.
- For batch queries, wrap this function in a loop and handle API rate limits.


#### **B. Retrieving Canonical Ensembl Transcript for a Mouse Gene**

This script queries the Ensembl REST API to:

- Look up a mouse gene by symbol.
- Retrieve the Ensembl canonical transcript ID as determined by Ensembl's internal prioritization.

```python
import requests

gene_symbol = "Kras"  # Mouse gene symbol, e.g., "Kras"
server = "https://rest.ensembl.org"
ext = f"/lookup/symbol/mus_musculus/{gene_symbol}?expand=1"

headers = {"Content-Type": "application/json"}
r = requests.get(server + ext, headers=headers)

if not r.ok:
    r.raise_for_status()

data = r.json()

# The 'canonical_transcript' field contains the Ensembl transcript ID designated as canonical
canonical_transcript_id = data.get('canonical_transcript')
print(f"Canonical Ensembl transcript for {gene_symbol}: {canonical_transcript_id}")
```

**Notes:**

- The canonical transcript in Ensembl is selected based on biological relevance and transcript completeness, and may differ from NCBI's choice.
- The transcript ID returned (e.g., ENSMUST000000XXXXXX) can be used to fetch sequence or further annotation via other Ensembl endpoints.
- For human, change `mus_musculus` to `homo_sapiens` in the URL.


#### **C. Additional Tips for Programmatic Access**

- **API Rate Limits:** Both NCBI and Ensembl APIs have request limits. For large-scale queries, add delays and monitor usage.
- **Annotation Versioning:** Always record the date and version of the database used for reproducibility.
- **Canonical Transcript Definitions:** The definition of "canonical" may differ between databases. For critical applications (e.g., clinical), consult database documentation or use gene-specific curation.
- **Batch Processing:** Both scripts can be adapted for batch processing by iterating over a list of gene symbols.

**References:**

- [NCBI E-utilities documentation](https://www.ncbi.nlm.nih.gov/books/NBK25501/)
- [Ensembl REST API documentation](https://rest.ensembl.org/documentation)

These code snippets and guidelines empower users to directly access and utilize canonical transcript IDs for downstream analyses, ensuring consistency and reproducibility in genomics workflows.



### **9. Additional Key Points and Resources**

- **Species-specific databases** (e.g., Wormbase, Flybase) and AWS iGenomes exist, but may be outdated[^1].
- **Chromosome naming conventions** differ by platform (e.g., NCBI/Ensembl: 1,2,...,X,Y,MT; UCSC: chr1, chr2, ..., chrM)[^1].
- **Annotation files should use stable gene IDs** (not gene symbols) for consistency across versions[^10].
- **Quality control**: Use annotation files with a `gene_biotype` field for best results in automated pipelines[^10].
- **Download locations**: NCBI, Ensembl, GENCODE, and UCSC all provide FTP and browser-based access to assemblies and annotations[^1][^5][^6][^7].


### **10. References to Key Sources**

- [^1] Introduction_to_resources_of_reference_genomes_for_omics_1744303221.pdf (attached PDF)
- [^8] Reference genome - Wikipedia
- [^4] PubMed: T2T-CHM13 reference assembly
- [^5] NCBI Insights: GRCm39
- [^2] UCSC News: Human pangenome
- [^9] NCBI: Annotating Genomes with GFF3 or GTF files
- [^6] NCBI Assembly: GRCm39
- [^3] Human Pangenome Reference Consortium
- [^11] Wikipedia: General Feature Format
- [^10] nf-core RNA-seq documentation
- [^12] Entrez Gene at NCBI
- [^7] NCBI Mouse Genome Overview


## **Summary Table: Human and Mouse Reference Resources**

| Species | Assembly | Annotation Source | Key FTP/Access Points | Notes |
| :-- | :-- | :-- | :-- | :-- |
| Human | GRCh38/hg38 | GENCODE, Ensembl, NCBI RefSeq | GENCODE, Ensembl, NCBI, UCSC | T2T-CHM13 and pangenome are most complete |
| Mouse | GRCm39/mm39 | GENCODE, Ensembl, NCBI RefSeq | GENCODE, Ensembl, NCBI, UCSC | GRCm38/mm10 for legacy data |

## **Final Notes**

- Always match annotation version to reference genome build.
- Use the most recent, well-curated annotation for new projects.
- For reproducibility, document all versions and sources used.
- No annotation is perfect; cross-check critical findings across databases when possible.

This guide integrates the most up-to-date, peer-reviewed, and platform-specific information for reference genome and annotation selection in human and mouse genomics[^1][^4][^5][^2][^3][^10][^7].

<div style="text-align: center">⁂</div>

[^1]: Introduction_to_resources_of_reference_genomes_for_omics_1744303221.pdf

[^2]: https://news.ucsc.edu/2023/05/pangenome-draft/

[^3]: https://humanpangenome.org

[^4]: https://pubmed.ncbi.nlm.nih.gov/38464973/

[^5]: https://ncbiinsights.ncbi.nlm.nih.gov/2020/10/16/refseq-grcm39/

[^6]: https://www.ncbi.nlm.nih.gov/assembly/7358741

[^7]: https://www.ncbi.nlm.nih.gov/grc/mouse

[^8]: https://en.wikipedia.org/wiki/Reference_genome

[^9]: https://www.ncbi.nlm.nih.gov/genbank/genomes_gff

[^10]: https://nf-co.re/rnaseq/3.16.1/docs/usage

[^11]: https://en.wikipedia.org/wiki/General_feature_format

[^12]: https://pmc.ncbi.nlm.nih.gov/articles/PMC1761442/

[^13]: https://www.ncbi.nlm.nih.gov/refseq/

[^14]: https://www.informatics.jax.org

[^15]: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.27/

[^16]: https://help.galaxyproject.org/t/updating-the-mouse-genome-to-grcm39-mm39-on-galaxy/5365

[^17]: https://hgdownload.soe.ucsc.edu

[^18]: https://www.biostars.org/p/231059/

[^19]: https://www.ncbi.nlm.nih.gov/genome/annotation_euk/Mus_musculus/109/

