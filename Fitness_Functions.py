import re
from Bio.Data import CodonTable
from Bio.Restriction import Analysis
from Bio.SeqUtils import GC
from Sequence_Container import SequenceContainer
from Bio_Structures import RibosomeBindingSites, splice_acceptors, splice_donors


def eval_host(individual, ancestor_sequence):
    assert (individual is SequenceContainer)
    sequence = getattr(individual, "sequence")
    return sum(['host'] for a, b in zip(sequence, ancestor_sequence) if a != b)


def eval_restriction_sites(individual, restrict_sites):
    """
    TODO: Make it remove rest sites
    """
    assert (individual is SequenceContainer)
    sequence = getattr(individual, "sequence")
    # check unwanted restriction sites
    analysis = Analysis(restrictionbatch=restrict_sites, sequence=sequence)
    result = analysis.full()
    # score the sequence based on the number of restriction sites
    score = 0
    for enz, cuts in result.items():
        for cut in cuts:
            score += 1
    return score


def eval_start_sites(individual, ribosome_binding_sites=RibosomeBindingSites, table_name="Standard"):
    """
    TODO: Make it remove start sites
    """
    assert (individual is SequenceContainer)
    sequence = getattr(individual, "sequence")
    codon_table = CodonTable.unambiguous_dna_by_name[table_name]

    # find all start codon sites (xTG)
    start_codon_positions = [
        m.start()
        for start_codon in codon_table.start_codons
        for m in re.finditer(start_codon, str(sequence))
    ]

    # None found
    if not len(start_codon_positions):
        return 0

    # check each start site for RBS
    # 18 base pairs upstream of each xTG, ignore 3 bp closest to xTG
    _rbs_offset = 18
    rbs_positions = [
        pos - _rbs_offset for pos in start_codon_positions if pos >= _rbs_offset
    ]
    mutable_seq = sequence.tomutable()

    score = 0
    for rbs_start in rbs_positions:
        # ignore 3 bp closest to xTG
        rbs_stop = rbs_start + _rbs_offset - 3
        rbs_query_seq = str(mutable_seq[rbs_start:rbs_stop])

        # check each unwanted RBS in each potential fragment
        for rbs, site in ribosome_binding_sites.items():
            search = rbs_query_seq.find(site)

            count = 0  # counter to prevent infinite loop
            while search != -1 and count < 10:
                # mutate residues if site is found
                codon_pos = (search + rbs_start) // 3
                for ii in range(2):
                    codon_idx = slice((codon_pos + ii) * 3, (codon_pos + ii + 1) * 3)
                    score += 1

                # reset sequence and search again
                rbs_query_seq = str(mutable_seq[rbs_start: rbs_stop + 3])
                search = rbs_query_seq.find(site)
                count += 1
    return score


def eval_repeats(individual, window_size=9):
    assert(individual is SequenceContainer)
    sequence = getattr(individual, "sequence")
    # iterate across overlapping chunks of complete codons
    codon_window = window_size // 3
    mutable_seq = sequence.tomutable()

    score = 0

    current_cycle = 0  # prevent infinite loops (caused by poly-TRP or poly-MET)
    keep_looping = True
    # `keep_looping` if any mutation is made,
    # i.e. continue until both checks pass without mutations
    while keep_looping and (current_cycle < (codon_window * 10)):
        keep_looping = False
        # iterate by codon, but map back to sequence-based indices
        for i in range(len(mutable_seq) // 3):
            window = slice(
                i * 3,
                (i + codon_window) * 3
                if (i + codon_window) * 3 < len(mutable_seq)
                else len(mutable_seq),
                )
            # make each mutable codon immutable so it can be hashed later
            codons = [
                str(mutable_seq[window][i: i + 3])
                for i in range(0, len(mutable_seq[window]), 3)
            ]
            # check if all codons in the window are identical
            if len(set(codons)) == 1:
                score += 1
            # check if the segment is found in the full sequence
            non_overlapping_matches = re.findall(
                str(mutable_seq[window]), str(mutable_seq)
            )
            if len(non_overlapping_matches) > 1 and len(mutable_seq[window]) > 3:
                score += len(non_overlapping_matches)
        current_cycle += 1
    return score


def eval_homopolymers(individual, n_codons=2, homopolymer_threshold=4):
    assert(individual is SequenceContainer)
    sequence = getattr(individual, "sequence")
    mutable_seq = sequence.tomutable()

    score = 0

    # look at each (n_codons * 3)-mer
    keep_looping = True
    while keep_looping:
        for i in range(0, len(mutable_seq), 3):
            window = slice(
                i,
                i + (n_codons * 3)
                if i + (n_codons * 3) < len(mutable_seq)
                else len(mutable_seq),
            )

            seq = str(mutable_seq[window])
            nt_count = 0
            nt_count_largest = 0
            last_letter = 'M'
            letter_largest = 'M'
            for letter in seq:
                if letter is last_letter:
                    nt_count += 1
                if nt_count > nt_count_largest:
                    nt_count_largest = nt_count
                    letter_largest = letter
                else:
                    # TODO: make loop stop if no larger sequence possible
                    last_letter = letter
                    nt_count = 1

                if nt_count_largest <= homopolymer_threshold:
                    keep_looping = False
                    continue

                score += nt_count // n_codons
                keep_looping = True
        return score


def eval_splice_sites(individual):
    assert(individual is SequenceContainer)
    sequence = getattr(individual, "sequence")

    def _pass_back_matches(list_of_sites, curr_dna):
        dna = str(curr_dna)
        sites = set(m for expr in list_of_sites for m in re.finditer(expr, dna))
        try:
            sites.remove(None)
        except KeyError:
            pass
        # remove redundancy
        sites = set((site.span(), site[0]) for site in sites)
        codon_bounds = [
            (s[0][0] // 3, -(-s[0][1] // 3)) for s in sorted(sites, key=lambda x: x[0])
        ]
        return codon_bounds

    def _get_splice_sites(curr_dna):
        donor_sites = _pass_back_matches(splice_donors, curr_dna)
        acceptor_sites = _pass_back_matches(splice_acceptors, curr_dna)
        return set(donor_sites + acceptor_sites)

    return len(_get_splice_sites(sequence.tomutable))


def eval_gc_content(individual, gc):
    assert(individual is SequenceContainer)
    sequence = getattr(individual, "sequence")
    window_size = gc.window_size  # tuples are immutable
    # some windows may be expressed as function of the sequence length
    if isinstance(window_size, str) and window_size.startswith("x"):
        window_size = int(float(window_size[1:]) * len(sequence))

    # iterate across overlapping chunks of complete codons
    codon_window = window_size // 3 + 1
    mutable_seq = sequence.tomutable()
    score = 0
    # iterate by codon, but map back to sequence-based indices
    for i in range(len(mutable_seq) // 3):
        window = slice(
            i * 3,
            (i + codon_window) * 3
            if (i + codon_window) * 3 < len(mutable_seq)
            else len(mutable_seq),
            )

        gc_percent = GC(mutable_seq[window]) / 100

        if gc_percent > gc.high:
            score += gc_percent - gc.high
        elif gc_percent < gc.low:
            score += gc.low - gc_percent
    return score


def eval_hairpins(individual, stem_length=10):
    """
    TODO: Make it remove hairpins
    :param individual:
    :param stem_length:
    :return:
    """
    assert(individual is SequenceContainer)
    sequence = getattr(individual, "sequence")
    mutable_seq = sequence.tomutable()
    score = 0
    for i in range(0, len(mutable_seq), stem_length):
        stem_seq = mutable_seq[i: i + stem_length].toseq()
        # include wobble base pairing for G-[CT]
        hairpin_pattern = "".join(
            [nt if nt != "C" else "[CT]" for nt in stem_seq.reverse_complement()]
        )
        for hairpin in re.finditer(hairpin_pattern, str(mutable_seq)):
            score += 1
    return score