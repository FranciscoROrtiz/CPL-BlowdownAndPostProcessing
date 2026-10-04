def AdiabaticBlowdown(V1, V2, P1, gamma):
    P2 = P1*(V1/V2)**gamma
    return P2
