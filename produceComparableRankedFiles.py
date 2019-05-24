from scipy import stats

def numberizeLine(parts):
    for i in range(1, len(parts)):
        if(parts[i] != "NA"):
            parts[i] = float(parts[i])
    return parts

def readMedians(filename, validIDs, nameConversion):
    genes = list()
    with open(filename) as medians:
        header = medians.readline().strip().split(" ")
        for line in medians:
            parts = line.strip().split(" ")
            if(nameConversion is None):
                parts[0] = parts[0].split(".")[0]
            elif(parts[0] in nameConversion):
                parts[0] = nameConversion[parts[0]]

            if(parts[0] in validIDs):
                genes.append(numberizeLine(parts))
    return (header, genes)

def readMapping(filename, delim, secondIndex, headerLine):
    firstToSecondMap = dict()
    secondToFirstMap = dict()
    alreadyEntered = set()
    with open(filename) as mapFile:
        if(headerLine):
            mapFile.readline()
        for line in mapFile:
            parts = line.strip().split(delim)
            if(len(parts) > secondIndex):
                assert(len(parts) == secondIndex + 1)
                if(parts[0] in alreadyEntered or parts[secondIndex] in alreadyEntered):
                    if(parts[0] in firstToSecondMap):
                        secondToFirstMap.pop(firstToSecondMap[parts[0]])
                        firstToSecondMap.pop(parts[0])

                    if(parts[secondIndex] in secondToFirstMap):
                        firstToSecondMap.pop(secondToFirstMap[parts[secondIndex]])
                        secondToFirstMap.pop(parts[secondIndex])

                    alreadyEntered.add(parts[0])
                    alreadyEntered.add(parts[secondIndex])
                else:
                    firstToSecondMap[parts[0]] = parts[secondIndex]
                    secondToFirstMap[parts[secondIndex]] = parts[0]
                    alreadyEntered.add(parts[0])
                    alreadyEntered.add(parts[secondIndex])
    return (firstToSecondMap, secondToFirstMap)

(dogNameConversion, notUsed) = readMapping("/seq/vgb/barkbase/LookupTablesStringtie/lookup.csv", ",", 3, False)
(dogToHumanMap, humanToDogMap) = readMapping("dog_to_human_orthologs_ensembl_96.31.txt", "\t", 1, True)

(dogHeader, dogGenes) = readMedians("dog-count-medians-output.txt", dogToHumanMap, dogNameConversion)
(humanHeader, humanGenes) = readMedians("human-count-medians-output.txt", humanToDogMap, None)
assert(dogHeader == humanHeader)

dogGeneMap = dict()
for currGene in dogGenes:
    dogGeneMap[currGene[0]] = currGene
humanGeneMap = dict()
for currGene in humanGenes:
    humanGeneMap[currGene[0]] = currGene

dogSet = set()
for currGene in dogGenes:
    dogSet.add(currGene[0])
commonGeneList = set()
for currGene in humanGenes:
    if(humanToDogMap[currGene[0]] in dogSet):
        commonGeneList.add(currGene[0])

corrMatrix = list()
for i in range(0, (len(dogHeader) - 1)*2):
    corrMatrix.append([-1] * ((len(dogHeader) - 1)*2))
pValueMatrix = list()
for i in range(0, (len(dogHeader) - 1)*2):
    pValueMatrix.append([-1] * ((len(dogHeader) - 1)*2))

for firstDogTissueIndex in range(1, len(dogHeader)):
    for secondDogTissueIndex in range(1, len(dogHeader)):
        firstCounts = list()
        secondCounts = list()
        for currGene in dogGenes:
            if(currGene[firstDogTissueIndex] > 1 or currGene[secondDogTissueIndex] > 1):
                firstCounts.append(currGene[firstDogTissueIndex])
                secondCounts.append(currGene[secondDogTissueIndex])
        assert(corrMatrix[firstDogTissueIndex - 1][secondDogTissueIndex - 1] == -1)
        assert(pValueMatrix[firstDogTissueIndex - 1][secondDogTissueIndex - 1] == -1)
        (corrMatrix[firstDogTissueIndex - 1][secondDogTissueIndex - 1], pValueMatrix[firstDogTissueIndex - 1][secondDogTissueIndex - 1]) = stats.spearmanr(firstCounts, secondCounts)

for firstHumanTissueIndex in range(1, len(humanHeader)):
    for secondHumanTissueIndex in range(1, len(humanHeader)):
        if(humanGenes[0][firstHumanTissueIndex] != "NA" and humanGenes[1][secondHumanTissueIndex] != "NA"):
            firstCounts = list()
            secondCounts = list()
            for currGene in humanGenes:
                if(currGene[firstHumanTissueIndex] > 1 or currGene[secondHumanTissueIndex] > 1):
                    firstCounts.append(currGene[firstHumanTissueIndex])
                    secondCounts.append(currGene[secondHumanTissueIndex])
            assert(corrMatrix[firstHumanTissueIndex - 1 + len(dogHeader) - 1][secondHumanTissueIndex - 1 + len(dogHeader) - 1] == -1)
            assert(pValueMatrix[firstHumanTissueIndex - 1 + len(dogHeader) - 1][secondHumanTissueIndex - 1 + len(dogHeader) - 1] == -1)
            (corrMatrix[firstHumanTissueIndex - 1 + len(dogHeader) - 1][secondHumanTissueIndex - 1 + len(dogHeader) - 1], pValueMatrix[firstHumanTissueIndex - 1 + len(dogHeader) - 1][secondHumanTissueIndex - 1 + len(dogHeader) - 1]) = stats.spearmanr(firstCounts, secondCounts)
        

for dogTissueIndex in range(1, len(dogHeader)):
    for humanTissueIndex in range(1, len(dogHeader)):
        if(humanGenes[0][humanTissueIndex] != "NA"):
            humanCounts = list()
            dogCounts = list()
            for currGene in commonGeneList:
                if(dogGeneMap[humanToDogMap[currGene]][dogTissueIndex] > 1 or humanGeneMap[currGene][humanTissueIndex] > 1):
                    dogCounts.append(dogGeneMap[humanToDogMap[currGene]][dogTissueIndex])
                    humanCounts.append(humanGeneMap[currGene][humanTissueIndex])
            assert(corrMatrix[dogTissueIndex - 1][humanTissueIndex - 1 + len(dogHeader) - 1] == -1)
            assert(pValueMatrix[dogTissueIndex - 1][humanTissueIndex - 1 + len(dogHeader) - 1] == -1)
            (corrMatrix[dogTissueIndex - 1][humanTissueIndex - 1 + len(dogHeader) - 1], pValueMatrix[dogTissueIndex - 1][humanTissueIndex - 1 + len(dogHeader) - 1]) = stats.spearmanr(dogCounts, humanCounts)

for dogTissueIndex in range(1, len(dogHeader)):
    for humanTissueIndex in range(1, len(dogHeader)):
        assert(corrMatrix[humanTissueIndex - 1 + len(dogHeader) - 1][dogTissueIndex - 1] == -1)
        assert(pValueMatrix[humanTissueIndex - 1 + len(dogHeader) - 1][dogTissueIndex - 1] == -1)
        (corrMatrix[humanTissueIndex - 1 + len(dogHeader) - 1][dogTissueIndex - 1], pValueMatrix[humanTissueIndex - 1 + len(dogHeader) - 1][dogTissueIndex - 1]) = (corrMatrix[dogTissueIndex - 1][humanTissueIndex - 1 + len(dogHeader) - 1], pValueMatrix[dogTissueIndex - 1][humanTissueIndex - 1 + len(dogHeader) - 1])

def printMatrix(matrix, header):
    tissueHeader = "\tdog_" + header[1]
    for i in range(2, len(header)):
        tissueHeader += "\tdog_" + header[i]
    for i in range(1, len(header)):
        if(matrix[0][i - 1 + len(header) - 1] != -1):
            tissueHeader += "\thuman_" + header[i]
    print(tissueHeader)
    for i in range(0, len(matrix)):
        if(i < len(header) - 1):
            currLine = "dog_" + header[i+1]
        else:
            currLine = "human_" + header[i + 1 - len(header) + 1]

        shouldPrint = False
        for k in range(0, len(matrix[i])):
            if(matrix[i][k] != -1):
                currLine += "\t" + str(matrix[i][k])
                shouldPrint = True
        if(shouldPrint):
            print(currLine)

printMatrix(corrMatrix, dogHeader)
#printMatrix(pValueMatrix, dogHeader)

                
