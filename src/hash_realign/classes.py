
class Sequence:
    def __init__(self, sequence):
        self.bases = sequence
        self.segSequence = []

    def add_sequence(self, sequence):
        sequence = sequence.upper()
        self.bases += sequence

    def clear(self):
        self.bases = ""

    def length(self):
        return len(self.bases)

    def get_bases(self):
        return self.bases

    def get_reverse_complement_bases(self):
        inv_seq = ""
        for i in range(len(self.bases) - 1, -1, -1):
            bp = self.bases[i]
            inv_bp = ''
            if bp == 'A':
                inv_bp = 'T'
            elif bp == 'T':
                inv_bp = 'A'
            elif bp == 'C':
                inv_bp = 'G'
            elif bp == 'G':
                inv_bp = 'C'
            else:
                inv_bp = 'N'

            inv_seq += inv_bp

        return inv_seq


class Segment:

    def __init__(self, read_start, ref_start, length, forward, seg_id):
        self.read_start = read_start
        self.ref_start = ref_start
        self.length = length
        self.forward = forward
        self.seg_id = seg_id
        if forward:
            self.read_end = self.read_start + (self.length - 1)
        else:
            self.read_end = self.read_start - (self.length - 1)
        self.ref_end = self.ref_start + (length - 1)


    def to_string(self):
        return str(self.read_start) + '\t' + str(self.read_end) + '\t' + str(self.ref_start) + '\t' + str(self.ref_end) + '\t' + str(self.forward)

