import time

time_list = []
len_list = (1, 10, 100, 1000)

for i in len_list[:2]:
    start = time.time()
    smis = ["CN(C)C(=O)c1ccccc1C(=O)N(C)C",
        "O=C1NC(=O)C(c2ccccc2)(c2ccccc2)N1",
        "CCCCNCc1ccc(S(=O)(=O)c2ccc(S(N)(=O)=O)s2)cc1",
        "COc1c(Cl)cc(Cl)c(OC)c1O",
        "Cc1cc(=O)[nH]c(=O)[nH]1",
        "O=c1[nH]c(=O)n(Cc2ccc(Cl)cc2)cc1F",
        "COc1ncnc2nccnc12",
        "CCC(CC)(Cc1ccc(C)cc1)C(=O)NO",
        "COc1ccc2cc(C(C)C(=O)OC3COC4C(O)COC34)ccc2c1",
        "CC(=O)Nc1cc(N[S+2]([O-])([O-])C(F)(F)F)c(C)cc1C",
        "O=C(O)c1cccnc1"] * i


    # feature time
    t1 = time.time() - start
    start = time.time()

    ## Run SDTNBI model
    from network.sdtnbi import SDTNBI
    from network.utils import save_network, load_pretrained_network

    # model = SDTNBI()
    # model.train(network="./data/drugbank-network.csv",
    #     drugs="./data/drugbank-drug.csv", r=2, nbits=2048)
    # save_network(model, "data/small-sdtnbi.pkl.gz") # uncomma for saving model

    model = load_pretrained_network("chembl")
    # model = load_pretrained_network("drugbank")

    F = model.predict(chemicals=smis, target='Q99640')   # Q99640 PMYT1_HUMAN

    # prediction time
    t2 = time.time() - start

    time_list.append((t1, t2))


for l, t in zip(len_list, time_list):
    print("Predict %d molecules:\tFeature Time: %.4f\tPrediction Time%.4f" % (l*10, t[0], t[1]))