# Python-Variant-Calling
Python implementation of viral trimming + variant calling, similar to iVar

Usage of call_variant_fetch.py:
python3 call_variant_fetch.py AlignmentFile ReferenceFile VariantOutputFilename
[-q minimal quality score to count base (Default:20)] [-t minimal frequency threshold to call variant (Default:0.03)]
[-m minimal number of reads to call variant (Default:0)]

Usage of call_variant_pileup.py:
python3 call_variant_pileup.py AlignmentFile ReferenceFile VariantOutputFilename
[-q minimal quality score to count base (Default:20)]
[-t minimal frequency threshold to call variant (Default:0.03)]
[-m minimal number of reads to call variant (Default:0)]
[-Q minimal quality score for a base to be considered during pileup (Default: 0)]
