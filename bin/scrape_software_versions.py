#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

regexes = {
    'qbicsoftware/rnadeseq': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'RColorBrewer': ['v_rcolorbrewer.txt', r"(\S+)"],
    'reshape2': ['v_reshape2.txt', r"(\S+)"],
    'Genefilter': ['v_genefilter.txt', r"(\S+)"],
    'DESeq2': ['v_deseq2.txt', r"(\S+)"],
    'ggplot2': ['v_ggplot2.txt', r"(\S+)"],
    'plyr': ['v_plyr.txt', r"(\S+)"],
    'vsn': ['v_vsn.txt', r"(\S+)"],
    'gplots': ['v_gplots.txt', r"(\S+)"],
    'pheatmap': ['v_pheatmap.txt', r"(\S+)"],
    'optparse': ['v_optparse.txt', r"(\S+)"],
    'svglite': ['v_optparse.txt', r"(\S+)"]
}
results = OrderedDict()
results['qbicsoftware/rnadeseq'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    with open(v[0]) as x:
        versions = x.read()
        match = re.search(v[1], versions)
        if match:
            results[k] = "v{}".format(match.group(1))

# Remove software set to false in results
for k in results:
    if not results[k]:
        del(results[k])

# Dump to YAML
print ('''
id: 'software_versions'
section_name: 'qbicsoftware/rnadeseq Software Versions'
section_href: 'https://github.com/qbicsoftware/rnadeseq'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in results.items():
    print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k,v))
print ("    </dl>")

# Write out regexes as tsv file:
with open('software_versions.tsv', 'w') as f:
    for k,v in results.items():
        f.write("{}\t{}\n".format(k,v))
