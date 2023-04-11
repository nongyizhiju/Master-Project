# This script would like to create a new feature for the certain gene
# insert the certain gene sequence into the reference genome
# update the positions of all downstream features and verify that the updated sequences and annotations are correct.

import gffutils
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from pyfaidx import Fasta

# Define the gene symbols of interest
genes_of_interest = ['OPN1LW']

# Load the reference genome and annotation file
genome_file = 'Homo_sapiens.GRCh38.dna.chromosome.X.fa'
annotation_file = 'Homo_sapiens.GRCh38.109.chromosome.X.gff3'

# Read the reference genome using pyfaidx
genome = Fasta(genome_file)

# Create a GFF3 database from the annotation file
db_file = 'annotation_Chro_X.db'
if os.path.exists(db_file):
    os.remove(db_file)

gffutils.create_db(annotation_file, db_file, id_spec=['ID', 'Name'], force=True, merge_strategy='create_unique')

# Open the GFF3 database
db = gffutils.FeatureDB(db_file)

# Define the output files
new_genome_file = 'duplicated_genes.fa'
new_annotation_file = 'duplicated_genes.gff3'

# Define the insert position
insert_position = 154159532
# insert the duplicated gene in 500bp downstream of the OPN1LW gene 

# Initialize a dictionary to store updated sequences
updated_genome = {seq_record.name: seq_record for seq_record in genome}

# Initialize the list of updated features
updated_features = []

for gene in db.features_of_type('gene'):
    if 'Name' not in gene.attributes:
        continue
    gene_symbol = gene.attributes['Name'][0]
    if gene_symbol not in genes_of_interest:
        continue

    # Calculate the length of the gene and the offset to apply to the duplicated gene and its child features
    gene_length = gene.end - gene.start + 1
    offset = gene_length

    # Extract the gene sequence
    gene_sequence = genome[gene.seqid][gene.start - 1:gene.end]

    # Insert the duplicated gene sequence into the reference genome
    genome_seq = str(genome[gene.seqid])
    new_genome_seq = genome_seq[:insert_position - 1] + str(gene_sequence) + genome_seq[insert_position - 1:]
    updated_genome[gene.seqid] = SeqRecord(Seq(new_genome_seq), id=gene.seqid)

    # Update the start and end positions of all downstream features in the GFF3 file
    for feature in db.region(seqid=gene.seqid, start=insert_position, end=len(new_genome_seq), featuretype=('gene', 'exon', 'CDS')):
        updated_feature = gffutils.Feature(
            seqid=feature.seqid,
            source=feature.source,
            featuretype=feature.featuretype,
            start=feature.start + offset,
            end=feature.end + offset,
            strand=feature.strand,
            attributes=dict(feature.attributes)
        )
        updated_feature.id = feature.id
        updated_features.append(updated_feature)

    # Create a new feature for the duplicated gene
    duplicated_gene = gffutils.Feature(
    seqid=gene.seqid,
    source=gene.source,
    featuretype=gene.featuretype,
    start=insert_position,
    end=insert_position + gene_length - 1,
    strand=gene.strand,
    attributes=dict(gene.attributes)
)
duplicated_gene.attributes['Name'] = [f"{gene_symbol}_dup"]
 

# Write the new FASTA file with the inserted duplicated gene
with open(new_genome_file, 'w') as new_genome_handle:
    for seq_record in updated_genome.values():
        SeqIO.write(seq_record, new_genome_handle, 'fasta')

# Write the new annotation file with updated features
with open(new_annotation_file, 'w') as output_handle:
    # Iterate through all original features in the GFF3 database
    for feature in db.all_features():
        # If the feature has been updated, write the updated feature instead
        updated_feature = next((uf for uf in updated_features if uf.id == feature.id), None)
        if updated_feature:
            output_handle.write(str(updated_feature) + '\n')
        else:
            output_handle.write(str(feature) + '\n')


# Verification steps -- check the inserted sequence and modified annotation file 
# 1. Compare the lengths of the original and new genome sequences
original_length = len(str(genome[gene.seqid]))
new_length = len(str(updated_genome[gene.seqid].seq))
print(f"Original genome length: {original_length}")
print(f"New genome length: {new_length}")
assert new_length == original_length + gene_length 
print('The length of the inserted sequence:', new_length - original_length)

# 2. Extract the inserted gene sequence and compare it to the original gene sequence
inserted_gene_sequence = updated_genome[gene.seqid].seq[insert_position - 1:insert_position - 1 + gene_length]
print(f"Original gene sequence (first 100 bp): {gene_sequence[:100]}")
print(f"Original gene sequence (last 100 bp): {gene_sequence[-100:]}")
print(f"Inserted gene sequence (first 100 bp): {inserted_gene_sequence[:100]}")
print(f"Inserted gene sequence (last 100 bp): {inserted_gene_sequence[-100:]}")
assert str(inserted_gene_sequence) == str(gene_sequence)

# 3. check the annotation file 
# Read the original and new annotation files
with open(annotation_file, 'r') as original_annotation_handle:
    original_annotation = original_annotation_handle.readlines()

with open(new_annotation_file, 'r') as new_annotation_handle:
    new_annotation = new_annotation_handle.readlines()

# Create dictionaries to store features by gene name
original_features = {}
new_features = {}

# Helper function to parse a GFF3 line into a dictionary
def parse_gff3_line(line):
    fields = line.strip().split('\t')
    attributes = dict(attribute.split('=') for attribute in fields[8].split(';'))
    return {
        'seqid': fields[0],
        'source': fields[1],
        'featuretype': fields[2],
        'start': int(fields[3]),
        'end': int(fields[4]),
        'score': fields[5],
        'strand': fields[6],
        'phase': fields[7],
        'attributes': attributes
    }

# Parse the original annotation file
for line in original_annotation:
    if not line.startswith('#'):
        feature = parse_gff3_line(line)
        if 'Name' in feature['attributes']:
            gene_name = feature['attributes']['Name']
            if gene_name in genes_of_interest:
                if gene_name not in original_features:
                    original_features[gene_name] = []
                original_features[gene_name].append(feature)

# Parse the new annotation file
for line in new_annotation:
    if not line.startswith('#'):
        feature = parse_gff3_line(line)
        if 'Name' in feature['attributes']:
            gene_name = feature['attributes']['Name']
            if gene_name in genes_of_interest:
                if gene_name not in new_features:
                    new_features[gene_name] = []
                new_features[gene_name].append(feature)

# Compare the original and new features
for gene_name in genes_of_interest:
    print(f"Comparing features for gene {gene_name}")
    for i, (original_feature, new_feature) in enumerate(zip(original_features[gene_name], new_features[gene_name])):
        if original_feature["attributes"]["ID"] == new_feature["attributes"]["ID"]:
            continue  # Skip comparing the original feature with itself
        
        # Check the duplicated gene
        duplicated_feature = None
        for new_feature in new_features[gene_name + "_dup"]:
            if new_feature['start'] == insert_position:
                duplicated_feature = new_feature
                break

        if duplicated_feature:
            print(f"  Original feature {i}: {original_feature}")
            print(f"  Duplicated feature {i}: {duplicated_feature}")
            assert duplicated_feature['start'] == insert_position

        # Check the updated downstream gene
        for new_feature in new_features[gene_name]:
            if new_feature['start'] == original_feature['start'] + gene_length:
                updated_feature = new_feature
                break

        if updated_feature:
            print(f"  Original feature {i}: {original_feature}")
            print(f"  Updated feature {i}: {updated_feature}")
            assert updated_feature['start'] == original_feature['start'] + gene_length


