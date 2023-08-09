import time
from network.sdtnbi import SDTNBI
from network.utils import save_network, load_network, load_pretrained_network


if __name__ == "__main__":

    # model = SDTNBI()
    # model.train(network="./data/chembl-network.csv", drugs="./data/chembl-drug.csv")
    # save_network(model, "data/chembl-sdtnbi.pkl.gz")

    # model = SDTNBI()
    # model.train(network="./data/drugbank-network.csv", drugs="./data/drugbank-drug.csv")
    # save_network(model, "data/drugbank-sdtnbi.pkl.gz")

    model = SDTNBI()
    model.train(network="./data/all-target-network.csv", drugs="./data/all-target-drug.csv")
    save_network(model, "./data/all-target-sdtnbi.pkl.gz")

    # with open("./data/test.smi") as f:
    #     smis = f.read().split("\n")
    # if len(smis[-1]) < 2:
    #     smis.pop(-1)

    # # smis *= 100

    # model = load_pretrained_network("chembl")
    # F = model.predict(chemicals=smis, target='Q99835')

    # time.sleep(1)
    
    # model = load_pretrained_network("drugbank")
    # F = model.predict(chemicals=smis, target='Q99835')
    
    # # print(F)


    with open("./data/test.smi") as f:
        smis = f.read().split("\n")
    if len(smis[-1]) < 2:
        smis.pop(-1)

    # smis *= 100

    model = load_network("./data/all-target-sdtnbi.pkl.gz")
    F = model.predict(chemicals=smis, target='Q99835')
    
    print(F)