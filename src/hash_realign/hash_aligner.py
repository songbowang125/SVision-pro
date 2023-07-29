#!/usr/bin/env python3


from src.hash_realign.classes import Segment


class HashAligner():

    def __init__(self, k, windowSize, mismatchNum, repeat_thresh):
        self.k = k
        self.windowSize = windowSize
        self.mismatchNum = mismatchNum

        self.segments = []
        self.selfDiffSegs = []
        self.compareDiffSegs = []
        self.N = 'N'
        self.G = 'G'
        self.A = 'A'
        self.T = 'T'
        self.C = 'C'
        self.avoid_kmers = []
        self.repeat_thresh = repeat_thresh
    def getSegments(self):
        return self.segments

    def getSelfDiffSegs(self):
        return self.selfDiffSegs

    def run(self, x, y, compareDiffSegs=None, y_hashvalue=None, avoid_kmers_from_ref=None):
        self.ref_length = len(y.get_bases())
        self.compareDiffSegs = compareDiffSegs
        self.y_hashvalues = y_hashvalue

        self.makePairwiseAlignment(x, y, avoid_kmers_from_ref)

    def extendKmersForward(self, xBases, yBases, matchPositions, p, i, segId):
        matchLength = self.k
        mismatch = 0

        # print(start, end)
        while mismatch <= self.mismatchNum:

            # beyond the length of x
            if matchPositions[p] + matchLength >= len(xBases) - 1:
                break
            # beyond the length of y
            if i + matchLength >= len(yBases) - 1:
                break

            xBase = xBases[matchPositions[p] + matchLength]
            yBase = yBases[i + matchLength]

            # stop extension when meet 'N'
            if xBase == self.N or yBase == self.N:
                break

            # stop extension when different
            if xBase != yBase:
                mismatch += 1

            matchLength += 1

        # when longer than windowSize
        if matchLength >= self.windowSize:
            d = Segment(matchPositions[p], i, matchLength, True, segId)

            if self.compareDiffSegs == None:
                self.segments.append(d)
                # if seq is ref(compareDiff = None), then calculate the ref's diff list
                if self.calDiffForRef(d):
                    self.selfDiffSegs.append(d)
            else:
                # all seg is diff

                if self.compareWithDiffSegs(d) is False:
                    # if d.get_read_start == 0 and d.get_ref_start == 0:
                    #     print(d.to_string())
                    self.segments.append(d)

    def extendKmersReverse(self, reverseXbases, yBases, matchPosition, i, segId):
        matchLength = self.k
        mismatch = 0
        while mismatch <= self.mismatchNum:
            if matchPosition + matchLength >= len(reverseXbases) - 1:
                break
            if i + matchLength >= len(yBases) - 1:
                break
            xBase = reverseXbases[matchPosition + matchLength]
            yBase = yBases[i + matchLength]

            if xBase == self.N or yBase == self.N:
                break

            # stop extension when different
            if xBase != yBase:
                mismatch += 1

            matchLength += 1
        if matchLength >= self.windowSize:
            d = Segment((len(reverseXbases) - 1) - matchPosition, i, matchLength, False, segId)

            if self.compareDiffSegs == None:
                self.segments.append(d)

                # if seq is ref(compareDiff = None), then calculate the ref's diff list
                # if not (d.read_start == 0 and d.ref_start == 0):
                #     self.selfDiffSegs.append(d)

                if self.calDiffForRef(d):
                    self.selfDiffSegs.append(d)
            else:
                if self.compareWithDiffSegs(d) is False:
                    self.segments.append(d)
                # if self.calDiff(d) or (abs(d.ref_end - d.ref_start) < y.length * 0.1):
                #     if self.compareWithDiffSegs(d) is False:
                #         self.segments.append(d)
                #
                # else:
                #     self.segments.append(d)

    def calHash(self, kmer):

        # hashValue = ''.join([str(bp) for bp in kmer])
        hashValue = kmer
        # hashValue = 0
        # for i in range(len(kmer)):
        #     hashValue += (int(kmer[9 - i] * math.pow(self.k, i)))
        #
        # hashValue = 0
        # power = 1
        # for d in range(9, -1, -1):
        #     if kmer[d] == self.G:
        #         hashValue += 0 * power
        #     elif kmer[d] == self.A:
        #         hashValue += 1 * power
        #     elif kmer[d] == self.T:
        #         hashValue += 2 * power
        #     elif kmer[d] == self.C:
        #         hashValue += 3 * power
        #
        #     power = power * 4
        return hashValue

    def makePairwiseAlignment(self, x, y, avoid_kmers_from_ref):

        xBases = x.get_bases()
        reverseXbases = x.get_reverse_complement_bases()
        # store hash value, format: hashvalue: [pos1, po2, ....]
        hashedPositions = {}

        # for sequence x, make hash
        for i in range(0, len(xBases) - (self.k + 1)):
            kmer = xBases[i: i + self.k]
            hashValue = self.calHash(kmer)

            if hashValue not in hashedPositions.keys():
                hashedPositions[hashValue] = []
            hashedPositions[hashValue].append(i)


        # get reverse seq, and make hash table
        for i in range(0, len(reverseXbases) - (self.k + 1)):
            kmer = reverseXbases[i: i + self.k]
            hashValue = self.calHash(kmer)

            if hashValue not in hashedPositions.keys():
                hashedPositions[hashValue] = []
            hashedPositions[hashValue].append(-1 - i)

        yBases = y.get_bases()
        segId = 0

        if self.y_hashvalues == None:

            self.hashvalues = []

            for i in range(0, len(yBases) - (self.k + 1)):
                kmer = yBases[i: i + self.k]
                hashValue = self.calHash(kmer)

                self.hashvalues.append(hashValue)

                if hashValue in hashedPositions.keys():
                    matchPositions = hashedPositions[hashValue]
                    # for each possible hit position
                    if len(matchPositions) >= self.repeat_thresh:
                        self.avoid_kmers.append(hashValue)
                    else:
                        for p in range(0, len(matchPositions)):

                            # hit in the same strand
                            # Kmer extension
                            if matchPositions[p] >= 0:
                                # Already matched in previous Kmer
                                if matchPositions[p] > 0 and i > 0 and xBases[matchPositions[p] - 1] == yBases[i - 1]:
                                    continue

                                self.extendKmersForward(xBases, yBases, matchPositions, p, i, segId)
                                segId += 1
                            else:
                                matchPosition = -1 - matchPositions[p]

                                if matchPosition > 0 and i > 0 and reverseXbases[matchPosition - 1] == yBases[i - 1]:
                                    continue

                                self.extendKmersReverse(reverseXbases, yBases, matchPosition, i, segId)
                                segId += 1
        else:
            for i in range(len(self.y_hashvalues)):
                hashValue = self.y_hashvalues[i]

                if hashValue in hashedPositions.keys():

                    if hashValue not in avoid_kmers_from_ref:

                        matchPositions = hashedPositions[hashValue]

                        # for each possible hit position
                        for p in range(0, len(matchPositions)):

                            # hit in the same strand
                            # Kmer extension
                            if matchPositions[p] >= 0:
                                # Already matched in previous Kmer
                                if matchPositions[p] > 0 and i > 0 and xBases[matchPositions[p] - 1] == yBases[i - 1]:
                                    continue

                                self.extendKmersForward(xBases, yBases, matchPositions, p, i, segId)
                                segId += 1
                            else:
                                matchPosition = -1 - matchPositions[p]

                                if matchPosition > 0 and i > 0 and reverseXbases[matchPosition - 1] == yBases[i - 1]:
                                    continue

                                self.extendKmersReverse(reverseXbases, yBases, matchPosition, i, segId)
                                segId += 1


    def getMergeSegments(self):
        return self.segments
        # main_segs = []
        # for seg in self.segments:
        #     if (seg.ref_start - 0) <= 5 and seg.forward == True:
        #         main_segs.append(seg)
        #     elif abs(seg.ref_end - self.ref_length) <= 5 and seg.forward == True:
        #         main_segs.append(seg)
        #
        # for main_seg in main_segs:
        #     self.segments.remove(main_seg)

        curSegNum = 1

        while curSegNum < len(self.segments):
            flag=0
            curSeg = self.segments[curSegNum]
            for i in range(curSegNum):
                candiSeg = self.segments[i]

                if self.linearOrNot(candiSeg, curSeg):
                    # print(candiSeg.read_start, candiSeg.read_end, candiSeg.ref_start, candiSeg.ref_end, curSeg.read_end, curSeg.ref_end)
                    # print('merge')
                    # merge and update segments list
                    if curSeg.forward == True:
                        candiSeg.read_end = max(curSeg.read_end, candiSeg.read_end)
                    elif curSeg.forward == False:
                        candiSeg.read_end = min(curSeg.read_end, candiSeg.read_end)
                        
                    candiSeg.ref_end = max(curSeg.ref_end, candiSeg.ref_end)

                    candiSeg.lenght = abs(candiSeg.length) + abs(curSeg.read_end - candiSeg.read_end)

                    # remove curSeg
                    self.segments.remove(curSeg)
                    flag = 1
                    break

            if flag == 0:
                curSegNum += 1

        # max_length = 0
        # for seg in self.segments:
        #     if seg.length > max_length:
        #         max_seg = seg
        #         max_length = seg.length

        # self.segments.extend(self.main_segs)
        after_filte_segs = []
        for seg in self.segments:
            if (seg.ref_end - seg.ref_start) >= 20:
                # print(seg.ref_end - seg.ref_start)
                after_filte_segs.append(seg)
        self.segments = after_filte_segs
        return self.segments


    def linearOrNot(self, i, j):

        # merge condition1: different strand, pass directly
        if i.forward != j.forward:
            return False

        # if i.forward is not False:
        #     return

        # merge condition2: colinear, not then false
        DIFF = self.calDiffBetTow(i, j)
        # print(DIFF)

        if DIFF > 1.2 or DIFF < 0.8:
            return False

        # merge condition3: colinear, yes, and dis
        DIS_X = abs(i.read_end - j.read_start)
        DIS_Y = abs(i.ref_end - j.ref_start)
        maxDis = (i.length + j.length) * 1.5
        if DIS_X > maxDis and DIS_Y > maxDis:
            return False

        # merge condition4: merged line's k is not -1 or 1
        tmp = float(j.read_end - i.read_start)
        if tmp == 0:
            tmp = 0.0001
        k = float(j.ref_end - i.ref_start) / tmp

        if abs((abs(k) - 1)) > 0.2:
            return False

        return True


    def compareWithDiffSegs(self, i):

        refStart = i.ref_start
        refEnd = i.ref_end
        forward = i.forward


        for tmpSeg in self.compareDiffSegs:
            # if tmpSeg.forward != forward:
            #     continue

            startDis = abs(refStart - tmpSeg.ref_start)
            endDis = abs(refEnd - tmpSeg.ref_end)

            if (startDis <= 5 and refEnd <= tmpSeg.ref_end) \
                    or (endDis <= 5 and refStart >= tmpSeg.ref_start):
                return True

        return False

    def calDiffForRef(self, i):
        # diff about the start
        diff2 = float(i.read_end) / float(i.ref_end)

        # diff about the center point
        centryX = float(i.read_start + i.read_end) / 2.0
        centryY = float(i.ref_start + i.ref_end) / 2.0
        diff3 = centryX / centryY

        if diff2 != 1 or diff3 != 1:
            return True
        else:
            return False

    def calDiff(self, i):
        diff1 = float(i.read_start + 1) / float(i.ref_start + 1)
        diff2 = float(i.read_end) / float(i.ref_end)

        centryX = float(i.read_start + i.read_end) / 2.0
        centryY = float(i.ref_start + i.ref_end) / 2.0

        diff3 = centryX / centryY

        thresh1 = 0.95
        thresh2 = 1.05

        if i.forward is True and diff1 == 1.0 and diff2 == 1.0 and diff3 == 1.0:
            return True
        else:
            return False
        # if i.forward is False and diff3 <= thresh2 and diff3 >= thresh1:
        #     return True
        #
        # if (diff1 > thresh2 or diff1 < thresh1) or (diff2 > thresh2 or diff2 < thresh1):
        #     return True
        #
        # else:
        #     return False

    def calDiffBetTow(self, i, j):
        if abs(float(i.ref_start - j.ref_start)) == 0:
            return 5
        return abs(float(i.read_start - j.read_start)) / abs(float(i.ref_start - j.ref_start))

