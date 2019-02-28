import csv
import os
import os.path
import sys
from itertools import islice

def loadParticipants(filename):
    participants = []
    with open(filename, "rb") as csvfile:
        reader = csv.reader(csvfile)
        # Skip first row, since it's a header row
        for (participantId, row) in enumerate(islice(reader, 1, None), 1):
            participantName = row[0]
            group = row[2]
            dtds = [int(d) for d in row[3:]]
            blocks = [dtds[:10], dtds[10:20], dtds[20:30], dtds[30:40]]
            participant = {'id': participantId, 'name': participantName, 'group': group, 'blocks': blocks}
            participants.append(participant)
    return participants

def loadSeqs(filename):
    nBlocks = 4
    seqsPerBlock = 10
    blockSeqs = []
    with open(filename, "rb") as seqfile:
        it = iter(seqfile)
        for b in range(nBlocks):
            block = []
            for i in range(seqsPerBlock):
                block.append(it.next())
            blockSeqs.append(block)
    return blockSeqs

def loadExcluded(filename):
    excluded = []
    with open(filename, "rb") as exfile:
        for line in exfile:
            if '#' in line:
                line = line[:line.find('#')]
            line = line.strip()
            if line:
                excluded.append(line)
    print excluded
    return excluded

def writeFile(target, participants, sequences):
    target.write("%d\n" % len(sequences))
    for s in sequences:
        target.write(s)
    target.write("%d\n" % len(participants))
    target.write("(1,%d) x (1,12)\n" % len(participants))
    target.write("[ ")
    for participant in participants:
        writeParticipant(participant['id'], participant['blocks'][block], target)
        target.write("\n")
    target.write("]\n")

def writeParticipant(participantNumber, draws, target):
    target.write(" %d %d %s" % (participantNumber, participantNumber + 1000, " ".join(str(d) for d in draws)))

if __name__ == "__main__":
    # Input files
    ALL_BLOCKS_DTD_FILENAME = "real_data/all_blocks_dtd.csv"
    ALL_BLOCKS_SEQ_FILENAME = "real_data/all_blocks_seq.csv"
    EXCLUDED_FILENAME = "real_data/excluded"

    output_dir = sys.argv[1] + "/real_input_data"

    VALID_GROUPS = ['0', '1', '2']

    # Start by loading all data, before writing anything (to avoid partial
    # overwrites if loading fails).
    participants = loadParticipants(ALL_BLOCKS_DTD_FILENAME)
    blockSeqs = loadSeqs(ALL_BLOCKS_SEQ_FILENAME)
    excluded = loadExcluded(EXCLUDED_FILENAME)

    participants = [p for p in participants
                    if p['group'] in VALID_GROUPS and p['name'] not in excluded]

    for participant in participants:
        groupName = ({'0': 'ARMS', '1': 'FEP', '2': 'Control', '66': 'Ideal'}
                [participant['group']])
        blockMeans = [sum(block)/len(block) for block in participant['blocks']]
        print "%s, %s, %s" % (participant['id'], groupName, blockMeans)
    print ''.join('===========\n' + ''.join(block) for block in blockSeqs)

    controls = [p for p in participants if p['group'] == '2']
    unhealthy = [p for p in participants if p['group'] == '1' or p['group'] == '0']
    for block in range(4):
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        dirname = output_dir + "/block%d" % (block + 1)
        if not os.path.isdir(dirname):
            os.mkdir(dirname)
        sequences = blockSeqs[block]
        with open(os.path.join(dirname, "control_data"), "w") as controlFile:
            writeFile(controlFile, controls, sequences)

        with open(os.path.join(dirname, "unhealthy_data"), "w") as unhealthyFile:
            writeFile(unhealthyFile, unhealthy, sequences)

        with open(os.path.join(dirname, "all_data"), "w") as allFile:
            writeFile(allFile, participants, sequences)


