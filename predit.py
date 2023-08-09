import time
import argparse

if __name__ == "__main__":
    from network.utils import load_pretrained_network
    from network.sdtnbi import logger
    start = time.time()

    parser = argparse.ArgumentParser(description='SDTNBI: an integrated network and chemoinformatics tool for systematic prediction of drug–target interactions and drug repositioning')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-smi', type=str, help='molecule in smiles format')
    group.add_argument('-txt', type=str, help='molecule in smiles list text')
    parser.add_argument('-t', type=str, required=True, help='target Uniprot ID')
    parser.add_argument('-o', type=str, required=False, help='Output file')
    parser.add_argument('-db', type=str, choices=["chembl", "drugbank"], default="drugbank", help='database')
    parser.add_argument('-v', action='store_true', help='print more logs')
    
    args = parser.parse_args()

    print(args.smi, args.txt)
    if args.smi:
        smis = [args.smi,]
    else:
        with open(args.txt) as f:
            smis = [l.strip().split()[0] for l in f.readlines()]

    # feature time
    t1 = time.time() - start

    ## Run SDTNBI model

    model = load_pretrained_network(args.db)
    if args.v:
        model.logger.setLevel(level=50)
    else:
        model.logger.setLevel(level=50)

    F = model.predict(chemicals=smis, target=args.t)   # Q99640 PMYT1_HUMAN

    # prediction time
    t2 = time.time() - start
    print("Predict %d molecules:\tFeature Time: %.4f\tPrediction Time%.4f" % (len(smis), t1, t2))

    results = sorted([x for x in zip(smis, F)], key=lambda x: x[1], reverse=True)
    # results = [x for x in zip(smis, F)]

    if args.o:
        with open(args.o, "w") as f:
            for smi, score in results:
                f.write("%s %f\n" % (smi, score))
    else:
        for smi, score in results:
            print("%s %f\n" % (smi, score))

