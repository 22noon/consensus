const CONFIG = {
  basePath: window.location.pathname.split('/')[1],
  "igvOptions": {
    "reference": {
      "id": "custom",
      "name": "MN632612.1",
      "fastaURL": "consensus.fa",
      "indexURL": "consensus.fa.fai"
    },
    "locus": "MN632612.1:4572",
    "tracks": [
      {
        "name": "Variants",
        "type": "variant",
        "format": "vcf",
        "url": "calls.norm.vcf.gz",
        "indexURL": "calls.norm.vcf.gz.tbi"
      }
    ]
  }
};