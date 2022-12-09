# mkTest
McDonald-Kreitman test written in Ruby. This code is quite old-- probably written in around 2002 -- so be careful!

# usage
To use this script you'll need a ruby interpretter (installed by default on some systems) and coding sequence data in 
fasta format. The input data must be in-frame (i.e. the first base represents the first codon position). 
Please see the provided example dataset.

```
$ ruby mkTest.rb
mkTest.rb ingroup.fa outgroup.fa
	options:
		-p outgroup2.fa (polarized MK test)
```

I've also provided a bit of test data (the venerable Adh locus from Drosophila melanogaster) to make sure
stuff works. So calling

```
ruby mkTest.rb kreitmanAdh.fa mauritianaAdh.fa
aaFix	aaPoly	silFix	silPoly
6	1	8	8
```

tells us that there were 6 amino acid fixations (aka nonsynonymous fixations), 1 amino acid polymorphism, 8 silent fixations,
and 8 silent polymorphisms. To get an associated p-value one can use a Fisher Exact test or similar on the resulting 2x2 table.

## options
In addition to the vanilla MK test, this script will do a polarized MK test along the ingroup lineage. That can be 
specified using the `-p` flag and then by providing a separate, second outgroup file. Polarization is done using 
a parsimony reconstruction of the ancestral state at each codon. 

# postscript
There is something like an easter egg in this script-- the first 3k lines or so are actually a set of Ruby classes that are usable for a large
number of population genetics tasks, way beyond MK tests. One could take those and do lots of other stuff, if you're interested. Sure, it makes for a bloated script but also pretty convenient :) 
