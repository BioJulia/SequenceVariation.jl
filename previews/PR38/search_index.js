var documenterSearchIndex = {"docs":
[{"location":"compare/","page":"Comparing variations","title":"Comparing variations","text":"CurrentModule = SequenceVariation","category":"page"},{"location":"compare/#Comparing-variations-in-sequences","page":"Comparing variations","title":"Comparing variations in sequences","text":"","category":"section"},{"location":"compare/#Checking-for-variations-in-a-known-haplotype","page":"Comparing variations","title":"Checking for variations in a known haplotype","text":"","category":"section"},{"location":"compare/","page":"Comparing variations","title":"Comparing variations","text":"Looking for a known Variation within a Haplotype is efficiently accomplished using the in operator.","category":"page"},{"location":"compare/","page":"Comparing variations","title":"Comparing variations","text":"using SequenceVariation, BioAlignments, BioSequences\n\nbovine = dna\"GACCGGCTGCATTCGAGGCTGCCAGCAAGCAG\";\novine  = dna\"GACCGGCTGCATTCGAGGCTGTCAGCAAACAG\";\nhuman  = dna\"GACAGGCTGCATCAGAAGAGGCCATCAAGCAG\";\n\nbos_ovis_alignment =\n    PairwiseAlignment(AlignedSequence(ovine, Alignment(\"32M\", 1, 1)), bovine);\nbos_human_alignment =\n    PairwiseAlignment(AlignedSequence(human, Alignment(\"32M\", 1, 1)), bovine);\n\nbos_ovis_haplotype = Haplotype(bos_ovis_alignment)\nbos_human_haplotype = Haplotype(bos_human_alignment)","category":"page"},{"location":"compare/","page":"Comparing variations","title":"Comparing variations","text":"println(\"\\tOvis aires\\tHomo sapiens\")\nfor v in vcat(variations(bos_ovis_haplotype), variations(bos_human_haplotype))\n    is_sheep = v in bos_ovis_haplotype\n    is_human = v in bos_human_haplotype\n    println(\"$v\\t$is_sheep\\t\\t$is_human\")\nend","category":"page"},{"location":"compare/#Constructing-new-haplotypes-based-on-other-variations","page":"Comparing variations","title":"Constructing new haplotypes based on other variations","text":"","category":"section"},{"location":"compare/","page":"Comparing variations","title":"Comparing variations","text":"New haplotypes can be constructed using variations. This might be useful to pool variations found on different reads or to filter variations from a haplotype that aren't validated by another haplotype.","category":"page"},{"location":"compare/","page":"Comparing variations","title":"Comparing variations","text":"sheeple = vcat(variations(bos_ovis_haplotype), variations(bos_human_haplotype));\nHaplotype(bovine, sheeple)\nreconstruct!(bovine, ans)","category":"page"},{"location":"api/","page":"API Reference","title":"API Reference","text":"CurrentModule = SequenceVariation\nDocTestSetup = quote\n    using SequenceVariation\nend","category":"page"},{"location":"api/#API-Reference","page":"API Reference","title":"API Reference","text":"","category":"section"},{"location":"api/#Edits","page":"API Reference","title":"Edits","text":"","category":"section"},{"location":"api/","page":"API Reference","title":"API Reference","text":"Substitution\nDeletion\nInsertion","category":"page"},{"location":"api/#SequenceVariation.Substitution","page":"API Reference","title":"SequenceVariation.Substitution","text":"Substitution\n\nRepresents the presence of a T at a given position. The position is stored outside this struct.\n\n\n\n\n\n","category":"type"},{"location":"api/#SequenceVariation.Deletion","page":"API Reference","title":"SequenceVariation.Deletion","text":"Deletion\n\nRepresents the deletion of N symbols. The location of the deletion is stored outside this struct\n\n\n\n\n\n","category":"type"},{"location":"api/#SequenceVariation.Insertion","page":"API Reference","title":"SequenceVariation.Insertion","text":"Insertion{S <: BioSequence}\n\nRepresents the insertion of a S into a sequence. The location of the insertion is stored outside the struct.\n\n\n\n\n\n","category":"type"},{"location":"api/#Variants","page":"API Reference","title":"Variants","text":"","category":"section"},{"location":"api/","page":"API Reference","title":"API Reference","text":"Haplotype\nreference(::Haplotype)\nvariations\nreconstruct\ntranslate(::Haplotype{S,T}, ::PairwiseAlignment{S,S}) where {S,T}","category":"page"},{"location":"api/#SequenceVariation.Haplotype","page":"API Reference","title":"SequenceVariation.Haplotype","text":"Haplotype{S<:BioSequence,T<:BioSymbol}\n\nA set of variations within a given sequence that are all found together. Depending on the field, it might also be referred to as a \"genotype\" or \"strain.\"\n\nConstructors\n\nHaplotype(ref::S, edits::Vector{Edit{S,T}}) where {S<:BioSequence,T<:BioSymbol}\nHaplotype(ref::S, vars::Vector{Variation{S,T}}) where {S<:BioSequence,T<:BioSymbol}\nHaplotype(\n    aln::PairwiseAlignment{T,T}\n) where {T<:LongSequence{<:Union{BS.AminoAcidAlphabet,BS.NucleicAcidAlphabet}}}\n\nWhen constructing a Haplotype from a vector of Edits or Variations, the edits are applied sequentially from first to last position, therefore the vector must always be sorted by position. These edits are sorted automatically if constructing from an alignment.\n\n\n\n\n\n","category":"type"},{"location":"api/#SequenceVariation.reference-Tuple{Haplotype}","page":"API Reference","title":"SequenceVariation.reference","text":"reference(h::Haplotype)\n\nGets the reference sequence of h.\n\n\n\n\n\n","category":"method"},{"location":"api/#SequenceVariation.variations","page":"API Reference","title":"SequenceVariation.variations","text":"variations(h::Haplotype{S,T}) where {S,T}\n\nConverts the Edits of h into a vector of Variations.\n\n\n\n\n\n","category":"function"},{"location":"api/#SequenceVariation.reconstruct","page":"API Reference","title":"SequenceVariation.reconstruct","text":"reconstruct(h::Haplotype)\n\nApply the edits in h to the reference sequence of h and return the mutated sequence\n\n\n\n\n\n","category":"function"},{"location":"api/#SequenceVariation.translate-Union{Tuple{T}, Tuple{S}, Tuple{Haplotype{S, T}, BioAlignments.PairwiseAlignment{S, S}}} where {S, T}","page":"API Reference","title":"SequenceVariation.translate","text":"translate(hap::Haplotype{S,T}, aln::PairwiseAlignment{S,S}) where {S,T}\n\nConvert the variations in hap to a new reference sequence based upon aln. The alignment rules follow the conventions of translate(::Variation, PairwiseAlignment). Indels at the beginning or end may not be preserved. Returns a new Haplotype\n\n\n\n\n\n","category":"method"},{"location":"api/#Variations","page":"API Reference","title":"Variations","text":"","category":"section"},{"location":"api/","page":"API Reference","title":"API Reference","text":"Variation\nreference(::Variation)\nmutation\ntranslate(::Variation{S,T}, ::PairwiseAlignment{S,S}) where {S,T}\nrefbases\naltbases","category":"page"},{"location":"api/#SequenceVariation.Variation","page":"API Reference","title":"SequenceVariation.Variation","text":"Variation{S<:BioSequence,T<:BioSymbol}\n\nA single change to a biological sequence. A general wrapper that can represent a sequence-specific Substitution, Deletion or Insertion. Variation is more robust than Edit, due to inclusion of the reference sequence and built-in validation.\n\nConstructors\n\nVariation(ref::S, e::Edit{S,T}) where {S<:BioSequence,T<:BioSymbol}\nVariation(ref::S, edit::AbstractString) where {S<:BioSequence}\n\nGenerally speaking, the Edit constructor should be avoided to ensure corectness: use of variations(::Haplotype) is encouraged, instead.\n\nConstructing a Variation from an AbstractString will parse the from edit using the following syntax:\n\nSubstitution: \"<REFBASE><POS><ALTBASE>\", e.g. \"G16C\"\nDeletion: \"Δ<STARTPOS>-<ENDPOS>\", e.g. \"Δ1-2\"\nInsertion: \"<POS><ALTBASES>\", e.g. \"11T\"\n\n\n\n\n\n","category":"type"},{"location":"api/#SequenceVariation.reference-Tuple{Variation}","page":"API Reference","title":"SequenceVariation.reference","text":"reference(v::Variation)\n\nGets the reference sequence of v\n\n\n\n\n\n","category":"method"},{"location":"api/#SequenceVariation.mutation","page":"API Reference","title":"SequenceVariation.mutation","text":"mutation(v::Variation)\n\nGets the underlying Substitution, Insertion, or Deletion of v.\n\n\n\n\n\n","category":"function"},{"location":"api/#SequenceVariation.translate-Union{Tuple{T}, Tuple{S}, Tuple{Variation{S, T}, BioAlignments.PairwiseAlignment{S, S}}} where {S, T}","page":"API Reference","title":"SequenceVariation.translate","text":"translate(var::Variation{S,T}, aln::PairwiseAlignment{S,S}) where {S,T}\n\nConvert the difference in var to a new reference sequence based upon aln. aln is the alignment of the old reference (aln.b) and the new reference sequence (aln.seq). Returns the new Variation.\n\n\n\n\n\n","category":"method"},{"location":"api/#SequenceVariation.refbases","page":"API Reference","title":"SequenceVariation.refbases","text":"refbases(v::Variation)\n\nGet the reference bases of v. Note that for deletions, refbases also returns the base before the deletion in accordance with the REF field of the VCF v4 specification.\n\n\n\n\n\n","category":"function"},{"location":"api/#SequenceVariation.altbases","page":"API Reference","title":"SequenceVariation.altbases","text":"altbases(v::Variation)\n\nGet the alternate bases of v. Note that for insertions, altbases also returns the base before the insertion in accordance with the ALT field of the VCF v4 specification.\n\n\n\n\n\n","category":"function"},{"location":"api/#Private-API","page":"API Reference","title":"Private API","text":"","category":"section"},{"location":"api/#Edits-2","page":"API Reference","title":"Edits","text":"","category":"section"},{"location":"api/","page":"API Reference","title":"API Reference","text":"Edit\n_mutation\n_lendiff","category":"page"},{"location":"api/#SequenceVariation.Edit","page":"API Reference","title":"SequenceVariation.Edit","text":"Edit{S <: BioSequence, T <: BioSymbol}\n\nAn edit of either Substitution{T}, Insertion{S} or Deletion at a position. If deletion: Deletion of length L at ref pos pos:pos+L-1 If insertion: Insertion of length L b/w ref pos pos:pos+1\n\n\n\n\n\n","category":"type"},{"location":"api/#SequenceVariation._mutation","page":"API Reference","title":"SequenceVariation._mutation","text":"_mutation(e::Edit)\n\nReturns the underlying Substitution, Insertion, or Deletion of e.\n\n\n\n\n\n","category":"function"},{"location":"api/#SequenceVariation._lendiff","page":"API Reference","title":"SequenceVariation._lendiff","text":"_lendiff(edit::Edit)\n\nGets the number of bases that edit adds to the reference sequence\n\n\n\n\n\n","category":"function"},{"location":"api/#Variants-2","page":"API Reference","title":"Variants","text":"","category":"section"},{"location":"api/","page":"API Reference","title":"API Reference","text":"_edits\n_is_valid(::Haplotype)","category":"page"},{"location":"api/#SequenceVariation._edits","page":"API Reference","title":"SequenceVariation._edits","text":"_edits(h::Haplotype)\n\nGets the Edits that comprise h\n\n\n\n\n\n","category":"function"},{"location":"api/#SequenceVariation._is_valid-Tuple{Haplotype}","page":"API Reference","title":"SequenceVariation._is_valid","text":"_is_valid(h::Haplotype{S,T}) where {S,T}\n\nValidate h. h is invalid if any of its operations are out of bounds, or the same position is affected by multiple edits.\n\n\n\n\n\n","category":"method"},{"location":"api/#Variations-2","page":"API Reference","title":"Variations","text":"","category":"section"},{"location":"api/","page":"API Reference","title":"API Reference","text":"_edit\n_is_valid(::Variation)","category":"page"},{"location":"api/#SequenceVariation._edit","page":"API Reference","title":"SequenceVariation._edit","text":"_edit(v::Variation)\n\nGets the underlying Edit of v\n\n\n\n\n\n","category":"function"},{"location":"api/#SequenceVariation._is_valid-Tuple{Variation}","page":"API Reference","title":"SequenceVariation._is_valid","text":"_is_valid(v::Variation)\n\nValidate v. v is invalid if its opertation is out of bounds.\n\n\n\n\n\n","category":"method"},{"location":"haplotypes/","page":"Working with haplotypes","title":"Working with haplotypes","text":"CurrentModule = SequenceVariation","category":"page"},{"location":"haplotypes/#Working-with-haplotypes","page":"Working with haplotypes","title":"Working with haplotypes","text":"","category":"section"},{"location":"haplotypes/#Calling-variants","page":"Working with haplotypes","title":"Calling variants","text":"","category":"section"},{"location":"haplotypes/","page":"Working with haplotypes","title":"Working with haplotypes","text":"The first step in working with sequence variation is to identify (call) variations between two sequences. SequenceVariation can directly call variants using the Haplotype(::PairwiseAlignment) constructor of the Haplotype type.","category":"page"},{"location":"haplotypes/","page":"Working with haplotypes","title":"Working with haplotypes","text":"using SequenceVariation, BioAlignments, BioSequences\n\nbovine = dna\"GACCGGCTGCATTCGAGGCTGCCAGCAAGCAG\";\novine  = dna\"GACCGGCTGCATTCGAGGCTGTCAGCAAACAG\";\nhuman  = dna\"GACAGGCTGCATCAGAAGAGGCCATCAAGCAG\";\n\nbos_ovis_alignment =\n    PairwiseAlignment(AlignedSequence(ovine, Alignment(\"32M\", 1, 1)), bovine);\nbos_human_alignment =\n    PairwiseAlignment(AlignedSequence(human, Alignment(\"32M\", 1, 1)), bovine);\n\nbos_ovis_haplotype = Haplotype(bos_ovis_alignment)\nbos_human_haplotype = Haplotype(bos_human_alignment)","category":"page"},{"location":"haplotypes/#Sequence-reconstruction","page":"Working with haplotypes","title":"Sequence reconstruction","text":"","category":"section"},{"location":"haplotypes/","page":"Working with haplotypes","title":"Working with haplotypes","text":"If the alternate sequence of a haplotype is no longer available (as is often the case when calling variants from alignment files), then the sequence can be retrieved using the reconstruct function.","category":"page"},{"location":"haplotypes/","page":"Working with haplotypes","title":"Working with haplotypes","text":"human2 = reconstruct(bos_human_haplotype)\nhuman2 == bovine\nhuman2 == human","category":"page"},{"location":"haplotypes/#Reference-switching","page":"Working with haplotypes","title":"Reference switching","text":"","category":"section"},{"location":"haplotypes/","page":"Working with haplotypes","title":"Working with haplotypes","text":"All variations within a haplotype can be mapped to a new reference sequence given an alignment between the new and old references using the translate function. This could be useful if variants were called against a reference sequence for the entire species, but need to be analyzed as variants of a subtype later.","category":"page"},{"location":"haplotypes/","page":"Working with haplotypes","title":"Working with haplotypes","text":"ovis_human_alignment =\n    PairwiseAlignment(AlignedSequence(human, Alignment(\"32M\", 1, 1)), ovine)\nSequenceVariation.translate(bos_ovis_haplotype, ovis_human_alignment)","category":"page"},{"location":"#SequenceVariation.jl","page":"Home","title":"SequenceVariation.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"(Image: Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.) (Image: Latest Release) (Image: MIT license) (Image: Stable documentation) (Image: Latest documentation)","category":"page"},{"location":"","page":"Home","title":"Home","text":"This project follows the semver pro forma and uses the OneFlow branching model.","category":"page"},{"location":"#Description","page":"Home","title":"Description","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"SequenceVariation provides a Julia vocabulary for comparing genetic mutations within biological sequences.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"You can install SequenceVariation from the Julia REPL. Press ] to enter pkg mode, and enter the following:","category":"page"},{"location":"","page":"Home","title":"Home","text":"add SequenceVariation","category":"page"},{"location":"#Testing","page":"Home","title":"Testing","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"SequenceVariation is tested against Julia 1.6 (LTS), 1.X (release) and nightly on Linux.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Latest build status:","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: Unit Tests) (Image: Documentation) (Image: codecov) (Image: Aqua QA)","category":"page"},{"location":"#Contributing","page":"Home","title":"Contributing","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"We appreciate contributions from users including reporting bugs, fixing issues, improving performance and adding new features.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Take a look at the contributing files for detailed contributor and maintainer guidelines, and code of conduct.","category":"page"},{"location":"#Financial-contributions","page":"Home","title":"Financial contributions","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"We also welcome financial contributions in full transparency on our open collective. Anyone can file an expense. If the expense makes sense for the development the core contributors and the person who filed the expense will be reimbursed.","category":"page"},{"location":"#Backers-and-Sponsors","page":"Home","title":"Backers & Sponsors","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Thank you to all our backers and sponsors!","category":"page"},{"location":"","page":"Home","title":"Home","text":"Love our work and community? Become a backer.","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: backers)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Does your company use BioJulia? Help keep BioJulia feature rich and healthy by sponsoring the project. Your logo will show up here with a link to your website.","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: ) (Image: ) (Image: ) (Image: ) (Image: ) (Image: ) (Image: ) (Image: ) (Image: ) (Image: )","category":"page"},{"location":"#Questions?","page":"Home","title":"Questions?","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"If you have a question about contributing or using BioJulia software, come on over and chat to us on the Julia Slack workspace, or you can try the Bio category of the Julia discourse site.","category":"page"},{"location":"variations/","page":"Working with variations","title":"Working with variations","text":"CurrentModule = SequenceVariation","category":"page"},{"location":"variations/#Working-with-individual-variations","page":"Working with variations","title":"Working with individual variations","text":"","category":"section"},{"location":"variations/#Construction","page":"Working with variations","title":"Construction","text":"","category":"section"},{"location":"variations/","page":"Working with variations","title":"Working with variations","text":"Individual Variations can be made using a reference sequence and string syntax","category":"page"},{"location":"variations/","page":"Working with variations","title":"Working with variations","text":"Variation type Syntax Interpretation Example\nSubstitutions <REF><POS><ALT> <ALT> is substituted for <REF> in position <POS> \"G16C\"\nDeletions Δ<START>-<END> All bases (inclusive) between <START> and <END> are deleted. It is valid to have <START> equal <END>: that is a deletion of one base. \"Δ1-2\"\nInsertions <POS><ALT> <ALT> is inserted between positions <POS> and <POS>+1 \"11T\"","category":"page"},{"location":"variations/","page":"Working with variations","title":"Working with variations","text":"using BioSequences: @dna_str\nusing SequenceVariation\nbovine_ins = dna\"GACCGGCTGCATTCGAGGCTGCCAGCAAGCAG\"\nVariation(bovine_ins, \"C4A\")\nmutation(ans)\ntypeof(mutation(Variation(bovine_ins, \"Δ13-14\")))\nmutation(Variation(bovine_ins, \"25ACA\"))","category":"page"},{"location":"variations/#Extraction","page":"Working with variations","title":"Extraction","text":"","category":"section"},{"location":"variations/","page":"Working with variations","title":"Working with variations","text":"Sequence variations may also be extracted wholesale from a Haplotype using the variations function.","category":"page"},{"location":"variations/","page":"Working with variations","title":"Working with variations","text":"using SequenceVariation, BioAlignments, BioSequences\n\nbovine = dna\"GACCGGCTGCATTCGAGGCTGCCAGCAAGCAG\";\novine  = dna\"GACCGGCTGCATTCGAGGCTGTCAGCAAACAG\";\nhuman  = dna\"GACAGGCTGCATCAGAAGAGGCCATCAAGCAG\";\n\nbos_ovis_alignment =\n    PairwiseAlignment(AlignedSequence(ovine, Alignment(\"32M\", 1, 1)), bovine);\nbos_human_alignment =\n    PairwiseAlignment(AlignedSequence(human, Alignment(\"32M\", 1, 1)), bovine);\n\nbos_ovis_haplotype = Haplotype(bos_ovis_alignment)\nbos_human_haplotype = Haplotype(bos_human_alignment)","category":"page"},{"location":"variations/","page":"Working with variations","title":"Working with variations","text":"variations(bos_ovis_haplotype)\nvariations(bos_human_haplotype)","category":"page"},{"location":"variations/#Reference-switching","page":"Working with variations","title":"Reference switching","text":"","category":"section"},{"location":"variations/","page":"Working with variations","title":"Working with variations","text":"An individual variation can be mapped to a new reference sequence given an alignment between the new and old references using the translate function.","category":"page"},{"location":"variations/","page":"Working with variations","title":"Working with variations","text":"ovis_human_alignment =\n    PairwiseAlignment(AlignedSequence(human, Alignment(\"32M\", 1, 1)), ovine)\nhuman_variation = first(variations(bos_ovis_haplotype))\nreference(ans) == bovine\nSequenceVariation.translate(human_variation, ovis_human_alignment)\nreference(ans) == bovine","category":"page"}]
}
