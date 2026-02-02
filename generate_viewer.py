#!/usr/bin/env python3
"""
Interactive Variant Viewer Generator
Generates an HTML-based IGV.js viewer for variants with read filtering and highlighting
"""

import pysam
import json
from jinja2 import Environment, FileSystemLoader
from pathlib import Path


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
    
    # Setup Jinja2 environment
    env = Environment(loader=FileSystemLoader(template_dir))
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
    vcf_path = "calls.norm.vcf.gz"
    bam_path = "alignment.markdup.bam"
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
