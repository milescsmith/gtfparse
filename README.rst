|GitHub| |Build Status| |Coverage Status| |CodeFactor| |Black|

gtfparse
========

Parsing tools for GTF (gene transfer format) files.

Example usage
=============

Parsing all rows of a GTF file into a Pandas DataFrame
------------------------------------------------------

.. code:: python

   from gtfparse import read_gtf

   # returns GTF with essential columns such as "feature", "seqname", "start", "end"
   # alongside the names of any optional keys which appeared in the attribute column
   df = read_gtf("gene_annotations.gtf")

   # filter DataFrame to gene entries on chrY
   df_genes = df[df["feature"] == "gene"]
   df_genes_chrY = df_genes[df_genes["seqname"] == "Y"]

Getting gene FPKM values from a StringTie GTF file
--------------------------------------------------

.. code:: python

   from gtfparse import read_gtf

   df = read_gtf(
       "stringtie-output.gtf",
       column_converters={"FPKM": float})

   gene_fpkms = {
       gene_name: fpkm
       for (gene_name, fpkm, feature)
       in zip(df["gene_name"], df["FPKM"], df["feature"])
       if feature == "gene"
   }

.. |GitHub| image:: https://img.shields.io/github/license/milescsmith/gtfparse
.. |Build Status| image:: https://travis-ci.com/milescsmith/gtfparse.svg?branch=master
   :target: https://travis-ci.com/milescsmith/gtfparse
.. |Coverage Status| image:: https://coveralls.io/repos/github/milescsmith/gtfparse/badge.svg?branch=master
   :target: https://coveralls.io/github/milescsmith/gtfparse?branch=master
.. |CodeFactor| image:: https://www.codefactor.io/repository/github/milescsmith/gtfparse/badge
   :target: https://www.codefactor.io/repository/github/milescsmith/gtfparse
.. |Black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/python/black
