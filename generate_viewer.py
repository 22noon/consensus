#!/usr/bin/env python3
"""
Interactive Variant Viewer Generator
Generates an HTML-based IGV.js viewer for variants with read filtering and highlighting
"""

import pysam
import json
from jinja2 import Environment, FileSystemLoader
from pathlib import Path
import argparse

parser = argparse.ArgumentParser(description="Generate IGV.js viewer")
parser.add_argument("--genome-name", help="Genome name/id")
parser.add_argument("--reference-fasta", help="Reference FASTA path/URL")
parser.add_argument("--bam-path", default="alignment.markdup.bam", help="BAM file path/URL")
parser.add_argument("--vcf-path", default="calls.norm.vcf.gz", help="VCF file path/URL")
args = parser.parse_args()
bam_path = args.bam_path
vcf_path = args.vcf_path

def read_variants_from_vcf(vcf_path, sequences):
    """Extract variants from VCF file"""
    vcf = pysam.VariantFile(vcf_path)
    variants = []
    Seq_in_VCF = []
    
    for record in vcf:
        if record.chrom not in Seq_in_VCF:
            Seq_in_VCF.append(record.chrom)
        variants.append({
            'chrom': record.chrom,
            'pos': record.pos,
            'ref': record.ref,
            'alt': ','.join([str(a) for a in record.alts]) if record.alts else '',
            'info': f"DP={record.info.get('DP', 'N/A')}"
        })
    
    vcf.close()
    print(f"Found {len(variants)} variants")

    for seq in sequences:
        if seq['name'] not in Seq_in_VCF:
            variants.append({
                'chrom': seq['name'],
                'pos': 0,
                'ref': '',
                'alt': '',
                'info': ""
            })

    return variants

def extract_sequence_names(bam_path):
    """Extract sequence names from BAM file header"""
    bam = pysam.AlignmentFile(bam_path, "rb")
    sequences = [{'name': ref['SN'], 'length': ref['LN']} for ref in bam.header['SQ']]
    bam.close()
    return sequences


def generate_html(variants, sequences, output_path, template_dir='templates'):
    """Generate HTML using Jinja2 template"""
    
    # Get first variant for initial locus
    first_locus = f"{variants[0]['chrom']}:{variants[0]['pos']}" if variants else "D388-WT_OR813926:1"
    
    reference_info = {
            'id': "custom",
            'name': "",
            'fastaURL': "",
            'indexURL': ""
    }

    # Update reference info from CLI
    reference_info["name"] = args.genome_name
    reference_info["fastaURL"] = args.reference_fasta
    reference_info["indexURL"] = f"{args.reference_fasta}.fai"

    # Set BAM/VCF from CLI (these may override later defaults)
    igv_options = {
        'reference': reference_info,
        'locus': first_locus,
        'tracks': [
            # {
            #     'name': 'All Alignments',
            #     'type': 'alignment',
            #     'format': 'bam',
            #     'url': bam_path,
            #     'indexURL': f"{bam_path}.bai",
            #     'height': 300,
            #     'viewAsPairs': False,
            #     'showCoverage': True
            #     # Note: 'filter' is set dynamically in JS via getFilterSettings()
            # },
            {
                'name': 'Variants',
                'type': 'variant',
                'format': 'vcf',
                'url': vcf_path,
                'indexURL': f"{vcf_path}.tbi"
            }
        ]
    }

    # Build config
    config = {
        'igvOptions': igv_options  # Add igvOptions to config
    }
    
    template_data = {
        'config_json': json.dumps(config, indent=2),
    }

    env = Environment(loader=FileSystemLoader(template_dir))
    # Render config.js
    config_template = env.get_template('config.js.j2')
    config_content = config_template.render(**template_data)
    with open('static/js/config.js', 'w') as f:
        f.write(config_content)
    
    # Setup Jinja2 environment
    template = env.get_template('viewer.html.j2')
    
    # Render template
    html_content = template.render(
        variants=variants,
        first_locus=first_locus,
        variants_json=json.dumps(variants),
        sequences_json=json.dumps(sequences)
    )
    
    # Write output
    with open(output_path, 'w') as f:
        f.write(html_content)
    
    print(f"Generated: {output_path}")


def main():
    """Main execution"""
    # Configuration
    output_path = "interactive_variants.html"

    print("Extracting sequences from BAM...")
    sequences = extract_sequence_names(bam_path)
    print(f"Found {len(sequences)} sequences")
    
    print("Reading variants from VCF...")
    variants = read_variants_from_vcf(vcf_path,sequences)

    print("Generating HTML...")
    generate_html(variants, sequences, output_path)

    print("Done!")


if __name__ == "__main__":
    main()
