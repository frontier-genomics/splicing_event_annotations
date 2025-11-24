#!/bin/sh
set -e

echo "Downloading Refseq (ncbiRefSeq)"
mkdir -p splicing_event_annotator/reference
docker run -e MYSQL_ALLOW_EMPTY_PASSWORD=1 --rm mysql mysql -ugenome -hgenome-euro-mysql.soe.ucsc.edu --compression-algorithms zlib -AD hg19 -BNe "SELECT r.bin,
    r.name,
    REPLACE(r.chrom, 'chr', '') AS chrom,
    r.strand,
    r.txStart,
    r.txEnd,
    r.cdsStart,
    r.cdsEnd,
    r.exonCount,
    CONVERT(r.exonStarts using utf8) AS exonStarts,
    CONVERT(r.exonEnds using utf8) AS exonEnds,
    r.score,
    r.name2,
    r.cdsStartStat,
    r.cdsEndStat,
    CONVERT(r.exonFrames using utf8) AS exonFrames,
    CASE WHEN r.name IN (SELECT name FROM hg19.ncbiRefSeqSelect) THEN 'Y' ELSE 'N' END AS MANE
FROM hg19.ncbiRefSeqCurated r
WHERE r.chrom IN ('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM')
;" > splicing_event_annotator/reference/hg19/genes.refGene

echo "Downloading Refseq (refGene)"
docker run -e MYSQL_ALLOW_EMPTY_PASSWORD=1 --rm mysql mysql -ugenome -hgenome-euro-mysql.soe.ucsc.edu --compression-algorithms zlib -AD hg19 -BNe "WITH ncbi_names AS
(
    SELECT DISTINCT REGEXP_REPLACE(name, '\.[0-9]+$', '') AS name 
    FROM hg19.ncbiRefSeqCurated
), version_info AS (
    SELECT acc, version
    FROM hgFixed.gbCdnaInfo
    WHERE organism = 3218 -- Homo sapien
    AND type = 'mRNA'
    AND SUBSTR(acc, 1, 2) IN ('NM','NR','XM','XR')
)
SELECT r.bin,
    CONCAT(r.name, '.', g.version) AS name,
    REPLACE(r.chrom, 'chr', '') AS chrom,
    r.strand,
    r.txStart,
    r.txEnd,
    r.cdsStart,
    r.cdsEnd,
    r.exonCount,
    CONVERT(r.exonStarts using utf8) AS exonStarts,
    CONVERT(r.exonEnds using utf8) AS exonEnds,
    r.score,
    r.name2,
    r.cdsStartStat,
    r.cdsEndStat,
    CONVERT(r.exonFrames using utf8) AS exonFrames,
    CASE WHEN r.name IN (SELECT name FROM hg19.ncbiRefSeqSelect) THEN 'Y' ELSE 'N' END AS MANE
FROM hg19.refGene r, version_info g
WHERE r.chrom IN ('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM')
AND r.name = g.acc
AND r.name NOT IN (SELECT name FROM ncbi_names)
;" >> splicing_event_annotator/reference/hg19/genes.refGene

echo "Done"