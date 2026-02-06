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

def read_variants_from_vcf(vcf_path):
    """Extract variants from VCF file"""
    vcf = pysam.VariantFile(vcf_path)
    variants = []
    
    for record in vcf:
        variants.append({
            'chrom': record.chrom,
            'pos': record.pos,
            'ref': record.ref,
            'alt': ','.join([str(a) for a in record.alts]) if record.alts else '',
            'info': f"DP={record.info.get('DP', 'N/A')}"
        })
    
    vcf.close()
    return variants


def extract_read_groups(bam_path):
    """Extract read groups from BAM file and assign colors"""
    rg_colors = [
        'rgb(255, 100, 100)',  # Red
        'rgb(100, 150, 255)',  # Blue
        'rgb(100, 200, 100)',  # Green
        'rgb(255, 200, 100)',  # Orange
        'rgb(200, 100, 255)',  # Purple
        'rgb(100, 200, 200)',  # Cyan
        'rgb(255, 150, 200)',  # Pink
        'rgb(150, 255, 100)',  # Light green
    ]
    
    bam = pysam.AlignmentFile(bam_path, "rb")
    
    # Add custom read group first
    read_groups = [
        {
            'id': 'NON_SPANNING_MATE',
            'sample': '',
            'library': 'unknown',
            'platform': 'ILLUMINA',
            'color': 'rgb(255, 100, 100)'
        },
        {
            'id': 'UNPAIRED',
            'sample': '',
            'library': 'unknown',
            'platform': 'ILLUMINA',
            'color': 'rgb(255, 100, 100)'
        },
        {
            'id': 'IMPROPER',
            'sample': '',
            'library': 'unknown',
            'platform': 'ILLUMINA',
            'color': 'rgb(255, 100, 100)'
        },
        {
            'id': 'PROPER',
            'sample': '',
            'library': 'unknown',
            'platform': 'ILLUMINA',
            'color': 'rgb(255, 100, 100)'
        },
        {
            'id': 'OVERLAP',
            'sample': '',
            'library': 'unknown',
            'platform': 'ILLUMINA',
            'color': 'rgb(255, 100, 100)'
        }
    ]
    
    # Extract read groups from BAM header
    color_index = 1
    if 'RG' in bam.header:
        for rg in bam.header['RG']:
            rg_dict = {
                'id': rg.get('ID', 'unknown'),
                'sample': rg.get('SM', 'unknown'),
                'library': rg.get('LB', 'unknown'),
                'platform': rg.get('PL', 'unknown'),
                'color': rg_colors[color_index % len(rg_colors)]
            }
            read_groups.append(rg_dict)
            color_index += 1
    
    bam.close()
    return read_groups


def generate_html(variants, read_groups, output_path, template_dir='templates'):
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
        read_groups=read_groups,
        first_locus=first_locus,
        variants_json=json.dumps(variants),
        read_groups_json=json.dumps(read_groups)
    )
    
    # Write output
    with open(output_path, 'w') as f:
        f.write(html_content)
    
    print(f"Generated: {output_path}")


def main():
    """Main execution"""
    # Configuration
    output_path = "interactive_variants.html"
    
    print("Reading variants from VCF...")
    variants = read_variants_from_vcf(vcf_path)
    print(f"Found {len(variants)} variants")
    
    print("Extracting read groups from BAM...")
    read_groups = extract_read_groups(bam_path)
    print(f"Found {len(read_groups)} read groups")
    
    print("Generating HTML...")
    generate_html(variants, read_groups, output_path)
    
    print("Done!")


if __name__ == "__main__":
    main()
