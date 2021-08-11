var documenterSearchIndex = {"docs":
[{"location":"#BGEN.jl","page":"BGEN.jl","title":"BGEN.jl","text":"","category":"section"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"Routines for reading compressed storage of genotyped or imputed markers","category":"page"},{"location":"#The-BGEN-Format","page":"BGEN.jl","title":"The BGEN Format","text":"","category":"section"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"Genome-wide association studies (GWAS) data with imputed markers are often saved in the BGEN format or .bgen file.","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"Used in:","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"Wellcome Trust Case-Control Consortium 2\nthe MalariaGEN project\nthe ALSPAC study\nUK Biobank: for genome-wide imputed genotypes and phased haplotypes","category":"page"},{"location":"#Features","page":"BGEN.jl","title":"Features","text":"","category":"section"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"Can store both hard-calls and imputed data\nCan store both phased haplotypes and phased genotypes\nEfficient variable-precision bit reapresntations\nPer-variant compression rightarrow easy to index \nSupported compression method: zlib and Zstandard. \nIndex files are often provided as .bgen.bgi files, which are plain SQLite3 databases.","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"Time to list variant identifying information (genomic location, ID and alleles): 18,496 samples, 121,668 SNPs (image source: https://www.well.ox.ac.uk/~gav/bgenformat/images/bgencomparison.png) (Image: )","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"Plink 1 format (.bed/.bim/.fam) has the list of variants as a separate file (.bim), effectively zero time.","category":"page"},{"location":"#Structure","page":"BGEN.jl","title":"Structure","text":"","category":"section"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"A header block followed by a series of [variant data block - (compressed) genotype data block] pairs.","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"Header block \nnumber of variants and samples\ncompression method (none, zlib or zstandard)\nversion of layout \nOnly \"layout 2\" is discussed below. \"Layout 1\" is also supported.\nsample identifiers (optional)","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"Variant data block\nvariant id\ngenomic position (chromosome, bp coordinate)\nlist of alleles\nGenotype data block (often compressed)\nploidy of each sample (may vary sample-by-sample)\nif the genotype data are phased\nprecision (B, number of bits to represent probabilities)\nprobabilitiy data (e.g. an unsigned B-bit integer x represents the probability of (fracx2^B-1)","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"BGEN.jl provides tools for iterating over the variants and parsing genotype data efficiently. It has been optimized for UK Biobank's zlib-compressed, 8-bit byte-aligned, all-diploid, all-biallelic datafiles.","category":"page"},{"location":"#Installation","page":"BGEN.jl","title":"Installation","text":"","category":"section"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"This package requires Julia v1.0 or later, which can be obtained from https://julialang.org/downloads/ or by building Julia from the sources in the https://github.com/JuliaLang/julia repository.","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"The package can be installed by running the following code:","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"using Pkg\npkg\"add BGEN\"","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"In order to run the examples below, the Glob package is also needed. ","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"pkg\"add Glob\"","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"versioninfo()","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"Julia Version 1.6.2\nCommit 1b93d53fc4 (2021-07-14 15:36 UTC)\nPlatform Info:\n  OS: macOS (x86_64-apple-darwin18.7.0)\n  CPU: Intel(R) Core(TM) i7-7820HQ CPU @ 2.90GHz\n  WORD_SIZE: 64\n  LIBM: libopenlibm\n  LLVM: libLLVM-11.0.1 (ORCJIT, skylake)","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"using BGEN, Glob","category":"page"},{"location":"#Example-Data","page":"BGEN.jl","title":"Example Data","text":"","category":"section"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"The example datafiles are stored in /data directory of this repository. It can be accessed through the function BGEN.datadir(). These files come from the reference implementation for the BGEN format.","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"Glob.glob(\"*\", BGEN.datadir())","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"78-element Vector{String}:\n \"/Users/xyz/.julia/dev/BGEN/src/../data/complex.10bits.bgen\"\n \"/Users/xyz/.julia/dev/BGEN/src/../data/complex.11bits.bgen\"\n \"/Users/xyz/.julia/dev/BGEN/src/../data/complex.12bits.bgen\"\n \"/Users/xyz/.julia/dev/BGEN/src/../data/complex.13bits.bgen\"\n \"/Users/xyz/.julia/dev/BGEN/src/../data/complex.14bits.bgen\"\n \"/Users/xyz/.julia/dev/BGEN/src/../data/complex.15bits.bgen\"\n \"/Users/xyz/.julia/dev/BGEN/src/../data/complex.16bits.bgen\"\n \"/Users/xyz/.julia/dev/BGEN/src/../data/complex.17bits.bgen\"\n \"/Users/xyz/.julia/dev/BGEN/src/../data/complex.18bits.bgen\"\n \"/Users/xyz/.julia/dev/BGEN/src/../data/complex.19bits.bgen\"\n \"/Users/xyz/.julia/dev/BGEN/src/../data/complex.1bits.bgen\"\n \"/Users/xyz/.julia/dev/BGEN/src/../data/complex.20bits.bgen\"\n \"/Users/xyz/.julia/dev/BGEN/src/../data/complex.21bits.bgen\"\n ⋮\n \"/Users/xyz/.julia/dev/BGEN/src/../data/example.6bits.bgen\"\n \"/Users/xyz/.julia/dev/BGEN/src/../data/example.7bits.bgen\"\n \"/Users/xyz/.julia/dev/BGEN/src/../data/example.8bits.bgen\"\n \"/Users/xyz/.julia/dev/BGEN/src/../data/example.8bits.bgen.bgi\"\n \"/Users/xyz/.julia/dev/BGEN/src/../data/example.9bits.bgen\"\n \"/Users/xyz/.julia/dev/BGEN/src/../data/example.gen\"\n \"/Users/xyz/.julia/dev/BGEN/src/../data/example.sample\"\n \"/Users/xyz/.julia/dev/BGEN/src/../data/example.v11.bgen\"\n \"/Users/xyz/.julia/dev/BGEN/src/../data/examples.16bits.bgen\"\n \"/Users/xyz/.julia/dev/BGEN/src/../data/haplotypes.bgen\"\n \"/Users/xyz/.julia/dev/BGEN/src/../data/haplotypes.bgen.bgi\"\n \"/Users/xyz/.julia/dev/BGEN/src/../data/haplotypes.haps\"","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"There are three different datasets with different format versions, compressions, or number of bits to represent probability values. ","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"example.*.bgen: imputed genotypes. \nhaplotypes.bgen: phased haplotypes. \ncomplex.*.bgen: includes imputed genotypes and phased haplotypes, and multiallelic genotypes.","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"Some of the .bgen files are indexed with .bgen.bgi files:","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"Glob.glob(\"*.bgen.bgi\", BGEN.datadir())","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"4-element Vector{String}:\n \"/Users/xyz/.julia/dev/BGEN/src/../data/complex.bgen.bgi\"\n \"/Users/xyz/.julia/dev/BGEN/src/../data/example.16bits.bgen.bgi\"\n \"/Users/xyz/.julia/dev/BGEN/src/../data/example.8bits.bgen.bgi\"\n \"/Users/xyz/.julia/dev/BGEN/src/../data/haplotypes.bgen.bgi\"","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"Sample identifiers may be either contained in the .bgen file, or is listed in an external .sample file.","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"Glob.glob(\"*.sample\", BGEN.datadir())","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"2-element Vector{String}:\n \"/Users/xyz/.julia/dev/BGEN/src/../data/complex.sample\"\n \"/Users/xyz/.julia/dev/BGEN/src/../data/example.sample\"","category":"page"},{"location":"#Type-Bgen","page":"BGEN.jl","title":"Type Bgen","text":"","category":"section"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"The type Bgen is the fundamental type for .bgen-formatted files. It can be created using the following line.","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"b = Bgen(BGEN.datadir(\"example.8bits.bgen\"); \n    sample_path=BGEN.datadir(\"example.sample\"), \n    idx_path=BGEN.datadir(\"example.8bits.bgen.bgi\"))","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"Bgen(IOStream(<file /Users/xyz/.julia/dev/BGEN/src/../data/example.8bits.bgen>), 0x000000000001f6ea, BGEN.Header(0x0000178c, 0x00000014, 0x000000c7, 0x000001f4, 0x01, 0x02, true), [\"sample_001\", \"sample_002\", \"sample_003\", \"sample_004\", \"sample_005\", \"sample_006\", \"sample_007\", \"sample_008\", \"sample_009\", \"sample_010\"  …  \"sample_491\", \"sample_492\", \"sample_493\", \"sample_494\", \"sample_495\", \"sample_496\", \"sample_497\", \"sample_498\", \"sample_499\", \"sample_500\"], Index(\"/Users/xyz/.julia/dev/BGEN/src/../data/example.8bits.bgen.bgi\", SQLite.DB(\"/Users/xyz/.julia/dev/BGEN/src/../data/example.8bits.bgen.bgi\"), UInt64[], String[], String[], UInt32[]))","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"The first argument is the path to the .bgen file. The optional keyword argument sample_path defines the location of  the .sample file. The second optional keyword argument idx_path determines the location of .bgen.bgi file. ","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"When a Bgen object is created, information in the header is parsed, and the index files are loaded if provided. You may retrieve basic information as follows. Variants are not yet parsed, and will be discussed later. ","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"io(b::Bgen): IOStream for the bgen file. You may also close this stream using close(b::Bgen).\nfsize(b::Bgen): the size of the bgen file.\nsamples(b::Bgen): the list of sample names. \nn_samples(b::Bgen): number of samples in the file.\nn_variants(b::Bgen): number of variants\ncompression(b::Bgen): the method each genotype block is compressed. It is either \"None\", \"Zlib\", or \"Zstd\".  ","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"io(b)","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"IOStream(<file /Users/xyz/.julia/dev/BGEN/src/../data/example.8bits.bgen>)","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"fsize(b)","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"128746","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"samples(b)","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"500-element Vector{String}:\n \"sample_001\"\n \"sample_002\"\n \"sample_003\"\n \"sample_004\"\n \"sample_005\"\n \"sample_006\"\n \"sample_007\"\n \"sample_008\"\n \"sample_009\"\n \"sample_010\"\n \"sample_011\"\n \"sample_012\"\n \"sample_013\"\n ⋮\n \"sample_489\"\n \"sample_490\"\n \"sample_491\"\n \"sample_492\"\n \"sample_493\"\n \"sample_494\"\n \"sample_495\"\n \"sample_496\"\n \"sample_497\"\n \"sample_498\"\n \"sample_499\"\n \"sample_500\"","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"n_samples(b)","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"500","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"n_variants(b)","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"199","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"compression(b)","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"\"Zlib\"","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"One may also access the list of RSIDs, chromosomes, and positions in chromosome of each variant stored using functions rsids(), chroms(), and positions(), respectively. ","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"rsids(b)","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"199-element Vector{String}:\n \"RSID_101\"\n \"RSID_2\"\n \"RSID_102\"\n \"RSID_3\"\n \"RSID_103\"\n \"RSID_4\"\n \"RSID_104\"\n \"RSID_5\"\n \"RSID_105\"\n \"RSID_6\"\n \"RSID_106\"\n \"RSID_7\"\n \"RSID_107\"\n ⋮\n \"RSID_194\"\n \"RSID_95\"\n \"RSID_195\"\n \"RSID_96\"\n \"RSID_196\"\n \"RSID_97\"\n \"RSID_197\"\n \"RSID_98\"\n \"RSID_198\"\n \"RSID_99\"\n \"RSID_199\"\n \"RSID_200\"","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"chroms(b)","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"199-element Vector{String}:\n \"01\"\n \"01\"\n \"01\"\n \"01\"\n \"01\"\n \"01\"\n \"01\"\n \"01\"\n \"01\"\n \"01\"\n \"01\"\n \"01\"\n \"01\"\n ⋮\n \"01\"\n \"01\"\n \"01\"\n \"01\"\n \"01\"\n \"01\"\n \"01\"\n \"01\"\n \"01\"\n \"01\"\n \"01\"\n \"01\"","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"positions(b)","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"199-element Vector{Int64}:\n   1001\n   2000\n   2001\n   3000\n   3001\n   4000\n   4001\n   5000\n   5001\n   6000\n   6001\n   7000\n   7001\n      ⋮\n  94001\n  95000\n  95001\n  96000\n  96001\n  97000\n  97001\n  98000\n  98001\n  99000\n  99001\n 100001","category":"page"},{"location":"#Variant-and-VariantIterator","page":"BGEN.jl","title":"Variant and VariantIterator","text":"","category":"section"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"As noted earlier, genotype information of each variant is compressed separately in .bgen files. The offsets (starting points in bgen file) of the genotypes may or may not be indexed by an external .bgen.bgi file. Thus, two ways to iterate over variants is provided through the function iterator(b; offsets=nothing, from_bgen_start=false). ","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"If offsets is provided, or .bgen.bgi is provided and ","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"from_bgen_start is false, it returns a VariantIteratorFromOffsets, iterating over the list of offsets.","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"Otherwise, it returns a VariantIteratorFromStart, iterating from the start of bgen file to the end of it sequentially. ","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"VariantIteratorFromOffsets and VariantIteratorFromStart are the subtypes of VariantIterator. ","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"Each element of VariantIterator is a Variant, containing the information of variants. We have following utility functions to access its information.","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"n_samples(v::Variant)\nvarid(v::Variant)\nrsid(v::Variant)\nchrom(v::Variant)\npos(v::Variant)\nn_alleles(v::Variant): number of alleles.\nalleles(v::Variant): list of alleles.","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"Merely the basic information of a variant is parsed for creating a Variant object. Nothing is decompressed, and genotype probabilities are not yet parsed yet. Decompression happens lazily, and is delayed until when we try to compute genotype probabilites or minor allele dosages (to be discussed later).","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"Since .bgen.bgi file is provided, the following order is based on the index file, sorted by genomic location.","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"for v in iterator(b) # \n    println(rsid(v))\nend","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"RSID_101\nRSID_2\nRSID_102\nRSID_3\nRSID_103\nRSID_4\nRSID_104\nRSID_5\nRSID_105\nRSID_6\nRSID_106\nRSID_7\nRSID_107\nRSID_8\nRSID_108\nRSID_9\nRSID_109\nRSID_10\nRSID_100\nRSID_110\nRSID_11\nRSID_111\nRSID_12\nRSID_112\nRSID_13\nRSID_113\nRSID_14\nRSID_114\nRSID_15\nRSID_115\nRSID_16\nRSID_116\nRSID_17\nRSID_117\nRSID_18\nRSID_118\nRSID_19\nRSID_119\nRSID_20\nRSID_120\nRSID_21\nRSID_121\nRSID_22\nRSID_122\nRSID_23\nRSID_123\nRSID_24\nRSID_124\nRSID_25\nRSID_125\nRSID_26\nRSID_126\nRSID_27\nRSID_127\nRSID_28\nRSID_128\nRSID_29\nRSID_129\nRSID_30\nRSID_130\nRSID_31\nRSID_131\nRSID_32\nRSID_132\nRSID_33\nRSID_133\nRSID_34\nRSID_134\nRSID_35\nRSID_135\nRSID_36\nRSID_136\nRSID_37\nRSID_137\nRSID_38\nRSID_138\nRSID_39\nRSID_139\nRSID_40\nRSID_140\nRSID_41\nRSID_141\nRSID_42\nRSID_142\nRSID_43\nRSID_143\nRSID_44\nRSID_144\nRSID_45\nRSID_145\nRSID_46\nRSID_146\nRSID_47\nRSID_147\nRSID_48\nRSID_148\nRSID_49\nRSID_149\nRSID_50\nRSID_150\nRSID_51\nRSID_151\nRSID_52\nRSID_152\nRSID_53\nRSID_153\nRSID_54\nRSID_154\nRSID_55\nRSID_155\nRSID_56\nRSID_156\nRSID_57\nRSID_157\nRSID_58\nRSID_158\nRSID_59\nRSID_159\nRSID_60\nRSID_160\nRSID_61\nRSID_161\nRSID_62\nRSID_162\nRSID_63\nRSID_163\nRSID_64\nRSID_164\nRSID_65\nRSID_165\nRSID_66\nRSID_166\nRSID_67\nRSID_167\nRSID_68\nRSID_168\nRSID_69\nRSID_169\nRSID_70\nRSID_170\nRSID_71\nRSID_171\nRSID_72\nRSID_172\nRSID_73\nRSID_173\nRSID_74\nRSID_174\nRSID_75\nRSID_175\nRSID_76\nRSID_176\nRSID_77\nRSID_177\nRSID_78\nRSID_178\nRSID_79\nRSID_179\nRSID_80\nRSID_180\nRSID_81\nRSID_181\nRSID_82\nRSID_182\nRSID_83\nRSID_183\nRSID_84\nRSID_184\nRSID_85\nRSID_185\nRSID_86\nRSID_186\nRSID_87\nRSID_187\nRSID_88\nRSID_188\nRSID_89\nRSID_189\nRSID_90\nRSID_190\nRSID_91\nRSID_191\nRSID_92\nRSID_192\nRSID_93\nRSID_193\nRSID_94\nRSID_194\nRSID_95\nRSID_195\nRSID_96\nRSID_196\nRSID_97\nRSID_197\nRSID_98\nRSID_198\nRSID_99\nRSID_199\nRSID_200","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"Setting from_bgen_start=true forces the iterator to iterate in the order of appearence in the bgen file. This may be different from the order in the index file.","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"for v in iterator(b; from_bgen_start=true)\n    println(rsid(v))\nend","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"RSID_2\nRSID_3\nRSID_4\nRSID_5\nRSID_6\nRSID_7\nRSID_8\nRSID_9\nRSID_10\nRSID_11\nRSID_12\nRSID_13\nRSID_14\nRSID_15\nRSID_16\nRSID_17\nRSID_18\nRSID_19\nRSID_20\nRSID_21\nRSID_22\nRSID_23\nRSID_24\nRSID_25\nRSID_26\nRSID_27\nRSID_28\nRSID_29\nRSID_30\nRSID_31\nRSID_32\nRSID_33\nRSID_34\nRSID_35\nRSID_36\nRSID_37\nRSID_38\nRSID_39\nRSID_40\nRSID_41\nRSID_42\nRSID_43\nRSID_44\nRSID_45\nRSID_46\nRSID_47\nRSID_48\nRSID_49\nRSID_50\nRSID_51\nRSID_52\nRSID_53\nRSID_54\nRSID_55\nRSID_56\nRSID_57\nRSID_58\nRSID_59\nRSID_60\nRSID_61\nRSID_62\nRSID_63\nRSID_64\nRSID_65\nRSID_66\nRSID_67\nRSID_68\nRSID_69\nRSID_70\nRSID_71\nRSID_72\nRSID_73\nRSID_74\nRSID_75\nRSID_76\nRSID_77\nRSID_78\nRSID_79\nRSID_80\nRSID_81\nRSID_82\nRSID_83\nRSID_84\nRSID_85\nRSID_86\nRSID_87\nRSID_88\nRSID_89\nRSID_90\nRSID_91\nRSID_92\nRSID_93\nRSID_94\nRSID_95\nRSID_96\nRSID_97\nRSID_98\nRSID_99\nRSID_100\nRSID_101\nRSID_102\nRSID_103\nRSID_104\nRSID_105\nRSID_106\nRSID_107\nRSID_108\nRSID_109\nRSID_110\nRSID_111\nRSID_112\nRSID_113\nRSID_114\nRSID_115\nRSID_116\nRSID_117\nRSID_118\nRSID_119\nRSID_120\nRSID_121\nRSID_122\nRSID_123\nRSID_124\nRSID_125\nRSID_126\nRSID_127\nRSID_128\nRSID_129\nRSID_130\nRSID_131\nRSID_132\nRSID_133\nRSID_134\nRSID_135\nRSID_136\nRSID_137\nRSID_138\nRSID_139\nRSID_140\nRSID_141\nRSID_142\nRSID_143\nRSID_144\nRSID_145\nRSID_146\nRSID_147\nRSID_148\nRSID_149\nRSID_150\nRSID_151\nRSID_152\nRSID_153\nRSID_154\nRSID_155\nRSID_156\nRSID_157\nRSID_158\nRSID_159\nRSID_160\nRSID_161\nRSID_162\nRSID_163\nRSID_164\nRSID_165\nRSID_166\nRSID_167\nRSID_168\nRSID_169\nRSID_170\nRSID_171\nRSID_172\nRSID_173\nRSID_174\nRSID_175\nRSID_176\nRSID_177\nRSID_178\nRSID_179\nRSID_180\nRSID_181\nRSID_182\nRSID_183\nRSID_184\nRSID_185\nRSID_186\nRSID_187\nRSID_188\nRSID_189\nRSID_190\nRSID_191\nRSID_192\nRSID_193\nRSID_194\nRSID_195\nRSID_196\nRSID_197\nRSID_198\nRSID_199\nRSID_200","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"With the presence of .bgen.bgi index file, one may select variants on a certain region using the function select_region(b, chrom; start=nothing, stop=nothing).  ","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"The following shows that all 199 variants in the bgen file are located on chromosome 01.  ","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"length(select_region(b, \"01\"))","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"199","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"We can see that the first variant since position 5000 at chromosome 01 is \"RSID_5\":  ","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"first(select_region(b, \"01\"; start=5000))","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"Variant(0x0000000000001ef8, 0x0000000000001f21, 0x0000000000002169, 0x00000248, 0x000001f4, \"SNPID_5\", \"RSID_5\", \"01\", 0x00001388, 0x0002, [\"A\", \"G\"], nothing)","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"And that the number of variants in chr01:5000-50000 is 92. ","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"length(select_region(b, \"01\"; start=5000, stop=50000))","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"92","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"Finally, one may use the parse_variants() function to retrieve the variant information as a Vector{Variant}. This is equivalent to calling collect() on the corresponding VariantIterator. It takes the same arguments as iterator(). This keeps all the information of variants in-memory. If the size of bgen file is too large, you might want to avoid this.","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"variants = parse_variants(b; from_bgen_start=true)","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"199-element Vector{Variant}:\n Variant(0x0000000000001790, 0x00000000000017b9, 0x0000000000001a82, 0x000002c9, 0x000001f4, \"SNPID_2\", \"RSID_2\", \"01\", 0x000007d0, 0x0002, [\"A\", \"G\"], nothing)\n Variant(0x0000000000001a82, 0x0000000000001aab, 0x0000000000001ced, 0x00000242, 0x000001f4, \"SNPID_3\", \"RSID_3\", \"01\", 0x00000bb8, 0x0002, [\"A\", \"G\"], nothing)\n Variant(0x0000000000001ced, 0x0000000000001d16, 0x0000000000001ef8, 0x000001e2, 0x000001f4, \"SNPID_4\", \"RSID_4\", \"01\", 0x00000fa0, 0x0002, [\"A\", \"G\"], nothing)\n Variant(0x0000000000001ef8, 0x0000000000001f21, 0x0000000000002169, 0x00000248, 0x000001f4, \"SNPID_5\", \"RSID_5\", \"01\", 0x00001388, 0x0002, [\"A\", \"G\"], nothing)\n Variant(0x0000000000002169, 0x0000000000002192, 0x0000000000002389, 0x000001f7, 0x000001f4, \"SNPID_6\", \"RSID_6\", \"01\", 0x00001770, 0x0002, [\"A\", \"G\"], nothing)\n Variant(0x0000000000002389, 0x00000000000023b2, 0x00000000000025df, 0x0000022d, 0x000001f4, \"SNPID_7\", \"RSID_7\", \"01\", 0x00001b58, 0x0002, [\"A\", \"G\"], nothing)\n Variant(0x00000000000025df, 0x0000000000002608, 0x00000000000027a4, 0x0000019c, 0x000001f4, \"SNPID_8\", \"RSID_8\", \"01\", 0x00001f40, 0x0002, [\"A\", \"G\"], nothing)\n Variant(0x00000000000027a4, 0x00000000000027cd, 0x00000000000029de, 0x00000211, 0x000001f4, \"SNPID_9\", \"RSID_9\", \"01\", 0x00002328, 0x0002, [\"A\", \"G\"], nothing)\n Variant(0x00000000000029de, 0x0000000000002a09, 0x0000000000002c43, 0x0000023a, 0x000001f4, \"SNPID_10\", \"RSID_10\", \"01\", 0x00002710, 0x0002, [\"A\", \"G\"], nothing)\n Variant(0x0000000000002c43, 0x0000000000002c6e, 0x0000000000002e8a, 0x0000021c, 0x000001f4, \"SNPID_11\", \"RSID_11\", \"01\", 0x00002af8, 0x0002, [\"A\", \"G\"], nothing)\n Variant(0x0000000000002e8a, 0x0000000000002eb5, 0x00000000000030e0, 0x0000022b, 0x000001f4, \"SNPID_12\", \"RSID_12\", \"01\", 0x00002ee0, 0x0002, [\"A\", \"G\"], nothing)\n Variant(0x00000000000030e0, 0x000000000000310b, 0x0000000000003375, 0x0000026a, 0x000001f4, \"SNPID_13\", \"RSID_13\", \"01\", 0x000032c8, 0x0002, [\"A\", \"G\"], nothing)\n Variant(0x0000000000003375, 0x00000000000033a0, 0x00000000000035dd, 0x0000023d, 0x000001f4, \"SNPID_14\", \"RSID_14\", \"01\", 0x000036b0, 0x0002, [\"A\", \"G\"], nothing)\n ⋮\n Variant(0x000000000001d991, 0x000000000001d9be, 0x000000000001dc12, 0x00000254, 0x000001f4, \"SNPID_189\", \"RSID_189\", \"01\", 0x00015ba9, 0x0002, [\"A\", \"G\"], nothing)\n Variant(0x000000000001dc12, 0x000000000001dc3f, 0x000000000001ddf2, 0x000001b3, 0x000001f4, \"SNPID_190\", \"RSID_190\", \"01\", 0x00015f91, 0x0002, [\"A\", \"G\"], nothing)\n Variant(0x000000000001ddf2, 0x000000000001de1f, 0x000000000001e011, 0x000001f2, 0x000001f4, \"SNPID_191\", \"RSID_191\", \"01\", 0x00016379, 0x0002, [\"A\", \"G\"], nothing)\n Variant(0x000000000001e011, 0x000000000001e03e, 0x000000000001e214, 0x000001d6, 0x000001f4, \"SNPID_192\", \"RSID_192\", \"01\", 0x00016761, 0x0002, [\"A\", \"G\"], nothing)\n Variant(0x000000000001e214, 0x000000000001e241, 0x000000000001e407, 0x000001c6, 0x000001f4, \"SNPID_193\", \"RSID_193\", \"01\", 0x00016b49, 0x0002, [\"A\", \"G\"], nothing)\n Variant(0x000000000001e407, 0x000000000001e434, 0x000000000001e6c9, 0x00000295, 0x000001f4, \"SNPID_194\", \"RSID_194\", \"01\", 0x00016f31, 0x0002, [\"A\", \"G\"], nothing)\n Variant(0x000000000001e6c9, 0x000000000001e6f6, 0x000000000001e8e1, 0x000001eb, 0x000001f4, \"SNPID_195\", \"RSID_195\", \"01\", 0x00017319, 0x0002, [\"A\", \"G\"], nothing)\n Variant(0x000000000001e8e1, 0x000000000001e90e, 0x000000000001ec86, 0x00000378, 0x000001f4, \"SNPID_196\", \"RSID_196\", \"01\", 0x00017701, 0x0002, [\"A\", \"G\"], nothing)\n Variant(0x000000000001ec86, 0x000000000001ecb3, 0x000000000001ef8b, 0x000002d8, 0x000001f4, \"SNPID_197\", \"RSID_197\", \"01\", 0x00017ae9, 0x0002, [\"A\", \"G\"], nothing)\n Variant(0x000000000001ef8b, 0x000000000001efb8, 0x000000000001f183, 0x000001cb, 0x000001f4, \"SNPID_198\", \"RSID_198\", \"01\", 0x00017ed1, 0x0002, [\"A\", \"G\"], nothing)\n Variant(0x000000000001f183, 0x000000000001f1b0, 0x000000000001f3d4, 0x00000224, 0x000001f4, \"SNPID_199\", \"RSID_199\", \"01\", 0x000182b9, 0x0002, [\"A\", \"G\"], nothing)\n Variant(0x000000000001f3d4, 0x000000000001f401, 0x000000000001f6ea, 0x000002e9, 0x000001f4, \"SNPID_200\", \"RSID_200\", \"01\", 0x000186a1, 0x0002, [\"A\", \"G\"], nothing)","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"If the index file (.bgi) is provided, the users may search for certain RSID in a BGEN file. ","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"v = variant_by_rsid(b, \"RSID_10\")","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"Variant(0x00000000000029de, 0x0000000000002a09, 0x0000000000002c43, 0x0000023a, 0x000001f4, \"SNPID_10\", \"RSID_10\", \"01\", 0x00002710, 0x0002, [\"A\", \"G\"], nothing)","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"Also, the users may look for the n-th (1-based) variant with respect to genomic location.","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"v = variant_by_index(b, 4)","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"Variant(0x0000000000001a82, 0x0000000000001aab, 0x0000000000001ced, 0x00000242, 0x000001f4, \"SNPID_3\", \"RSID_3\", \"01\", 0x00000bb8, 0x0002, [\"A\", \"G\"], nothing)","category":"page"},{"location":"#Genotype/haplotype-probabilities-and-minor-allele-dosage","page":"BGEN.jl","title":"Genotype/haplotype probabilities and minor allele dosage","text":"","category":"section"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"The genotype information is decompressed and parsed when probability data is needed. The parsing is triggered by a call to one of:","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"probabilities!(b::Bgen, v::Variant; T=Float64) : probability of each genotype/haplotype.\nminor_allele_dosage!(b::Bgen, v::Variant; T=Float64) : minor allele dosage for a biallelic variant.","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"Once parsed, the results are cached and loaded on any subsequent calls.  After that, one may access genotype information using the following functions, as well as probabilities!() and minor_allele_dosage!():","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"phased(v::Variant): if the stored data is phased\nmin_ploidy(v::Variant): minimum ploidy across the samples\nmax_ploidy(v::Variant): maximum ploidy across the samples\nploidy(v::Variant) : Vector of ploidy for each sample\nbit_depth(v::Variant) : number of bits used to represent a probability value\nmissings(v::Variant) : list of samples data is missing","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"These functions are allowed after calling minor_allele_dosage!():","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"minor_allele(v::Variant)\nmajor_allele(v::Variant)","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"If the data are not phased, probabilities!(b, v)[i, j] represents the probability of genotype i for sample j. Each column sums up to one. The genotypes are in colex-order of allele counts. For example, for three alleles with ploidy 3:","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"row index allele counts genotype\n1 (3, 0, 0) 111\n2 (2, 1, 0) 112\n3 (1, 2, 0) 122\n4 (0, 3, 0) 222\n5 (2, 0, 1) 113\n6 (1, 1, 1) 123\n7 (0, 2, 1) 223\n8 (1, 0, 2) 133\n9 (0, 1, 2) 233\n10 (0, 0, 3) 333","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"probabilities!(b, variants[1])","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"3×500 Matrix{Float32}:\n NaN  0.027451    0.0156863  0.0235294  …  0.0156863  0.921569   0.00392157\n NaN  0.00784314  0.0509804  0.933333      0.027451   0.0509804  0.984314\n NaN  0.964706    0.933333   0.0431373     0.956863   0.027451   0.0117647","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"Genotype data for sample 1 is missing in this case. ","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"missings(variants[1])","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"1-element Vector{Int64}:\n 1","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"On the other hand, if the data are phased, probabilities!(b, v)[i, j] represents the probability that haplotype (i - 1) ÷ n_alleles + 1 has allele (i - 1) % n_alleles + 1 for sample j, where n_alleles is the number of alleles. The below is an example of phased probabilities. i.e., each column represents each sample, and each group of n_alleles rows represent the allele probabilities for each haplotype. In this case, ploidy is [1, 2, 2, 2], thus indexes [3:4, 1] are invalid, and is filled with NaN.","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"b2 = Bgen(BGEN.datadir(\"complex.bgen\"))\nvs = parse_variants(b2)\np = probabilities!(b2, vs[3])","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"4×4 Matrix{Float32}:\n   1.0  0.0  1.0  1.0\n   0.0  1.0  0.0  0.0\n NaN    0.0  1.0  0.0\n NaN    1.0  0.0  1.0","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"This variant has two possible alleles (allele 1: \"A\" and allele 2: \"G\"), and all the samples are diploids except for the first one, which is monoploid.","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"It corresponds to a line of VCF file:","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  sample_0        sample_1        sample_2        sample_3\n01      3       V3      A       G       .       .       .       GT:HP   0:1,0   1|1:0,1,0,1     0|0:1,0,1,0     0|1:1,0,0,1","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"So the first sample is monoploid of A, the second sample is homozygote A|A, the third sample is homozygote G|G, and the last sample is heterozygote A|G (phased). ","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"We can confirm the phasedness and ploidy of each sample as follows.","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"phased(vs[3])","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"0x01","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"ploidy(vs[3])","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"4-element Vector{UInt8}:\n 0x01\n 0x02\n 0x02\n 0x02","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"alleles(vs[3])","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"2-element Vector{String}:\n \"A\"\n \"G\"","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"For biallelic nonphased genotype data, minor_allele_dosage!(b, v) can be computed. ","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"minor_allele_dosage!(b, variants[1])","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"500-element Vector{Float32}:\n NaN\n   0.0627451\n   0.08235294\n   0.9803922\n   0.09019608\n   0.14117648\n   1.0745099\n   0.054901965\n   0.10980393\n   0.121568635\n   0.14117648\n   0.21568629\n   0.08235294\n   ⋮\n   0.09411766\n   0.10196079\n   0.027450982\n   0.96470594\n   0.0\n   1.0117648\n   0.043137256\n   0.0627451\n   1.0431373\n   0.05882353\n   1.8941176\n   0.99215686","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"phased(variants[1])","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"0x00","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"n_alleles(variants[1])","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"2","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"minor_allele(variants[1])","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"\"A\"","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"major_allele(variants[1])","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"\"G\"","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"minor_allele_dosage!() supports a keyword argument mean_impute, which imputes missing value with the mean of the non-missing values.","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"minor_allele_dosage!(b, variants[1]; T=Float64, mean_impute=true)","category":"page"},{"location":"","page":"BGEN.jl","title":"BGEN.jl","text":"500-element Vector{Float32}:\n 0.39581063\n 0.0627451\n 0.08235294\n 0.9803922\n 0.09019608\n 0.14117648\n 1.0745099\n 0.054901965\n 0.10980393\n 0.121568635\n 0.14117648\n 0.21568629\n 0.08235294\n ⋮\n 0.09411766\n 0.10196079\n 0.027450982\n 0.96470594\n 0.0\n 1.0117648\n 0.043137256\n 0.0627451\n 1.0431373\n 0.05882353\n 1.8941176\n 0.99215686","category":"page"}]
}
