# Non Interactive Zero Knowledge Proofs of Graph 3- Colorings Writeup
#### By Jason Stillerman

## Problem Statement

Zero knowledge proofs are powerful tools that allow a Prover to convince a Verifier that the prover knows a witness, $w$, that satisfies $f(w) = True$ without revealing the witness itself. This turns out to be very useful for various application in cryptography, but on its own has a weakness: it requires many continuous interactions between the Prover and Verifier. This means that both parties have to be online during the verification process which is a limitation if you want to have many verifiers or communication is taking place on a high latency chanel. The number of rounds of communication is determined by how sure the Verifier needs to be that the Prover is not "getting lucky". In any given round of the proof, the Prover could be lying and getting away with it with some probability, $P$. It is only after $n$ rounds of this that you get $P^n$, the probability that Prover has gotten away with cheating after $n$ rounds, to be very low.

Adding non-interactivty to the ZKP allows the Prover to publish the proof in one step, and then after the fact as many verifiers as would like can verify the proof without interacting with the Prover. This is extremely useful, yet counter intuitive. The prover must pick their own challenges, which intuitively sounds untrustworthy, but there is a scheme for creating what is effectively a random oracle - a black box that spits out a random number the Prover cannot control, and we will discuss the technical details of that implementation below.

## The Naive Approach

The interactive ZKP protocol is as follows

1. Prover commits to a random coloring of G
2. Verifier picks an edge of G to challenge Prover on
3. Prover opens commitment of nodes on either side of the edge
4. Verifier confirms that colors are in fact different and commitments were respected
5. Repeat steps 1-4 until Verifier is satisfied that Prover did not just get lucky

To make this non-interactive, all of the Provers steps must happen before any of the Verifer(s)'s steps. To do this, we need a way for the Prover to pick the challenge edge for herself, but in a way that the Verifier can be convinced that it is actually random. To do this we will put the commitments through a good hash function and mod the result by the number of edges. Good hash functions are non-invertible and will map inputs uniformishly to the output space. The new protocol is as follows:

1. Prover commits to random coloring of G
2. Prover hashes that commitment and mods by the number of edges to get a random challenge edge
3. Prover reveals the commitments for that edge
4. Repeat steps 1-3 and publish the results each time
5. Verifier takes every example and
   1. Checks that the random oracle was operated fairly
   2. Checks that the commitments were respected
   3. Checks that the colors are distinct


Now the Prover can do all of its work offline, exchange only the bundle of commitments, and then the Verifier can do all of its verification offline as well. Interactivity removed!

## The Issue

Let's say the Prover does not actually have a valid 3-coloring of G; one of the edges has the same color on both sides of it. If the Prover picks a random recoloring, commits to it, hashes the commitment, gets a random challenge edge that doesn't reveal the cheating, they publish it. If it *does* reveal the cheating, they can pick a new random recoloring that hopefully does not end up choosing the cheating edge, and just not publish the commitment that revealed the cheating.

## The Solution

To combat this, we need some assurance that the Prover cannot pick and choose which commitments its going to publish. To achieve this, we make the Prover commit to all of the random colorings up front, and then generate a sequence of challenges based on *all* of the commitments. The library of commitments is hashed to get the first oracle, and that output is the input to the next oracle. We have effectively created a deterministic random number *sequence* oracle that can be verified.

If the Prover wants to cheat, they would need to find a huge set of commitments that when put into the random sequence oracle never reveals the challenge edge. This will be increasingly statistically unlikely for larger and larger numbers of commitments and a good hash function. The adapted protocol is as follows:

1. Prover commits to N random colorings of G
2. Prover hashes *all* commitments and gets the seed for the random oracle
3. Take the output of the oracle and feed it back into itself to get a new random number. Repeat this until you have N random numbers we will use as the challenge edges
4. Prover opens the commitment for each challenge edge for each commitment
5. Prover broadcasts all commitments, challenge edges, and open commitments.
6. Verifier takes every example and
   1. Checks that the random oracle was operated fairly
   2. Checks that the commitments were respected
   3. Checks that the colors are distinct


## Results
Assuming there is only a single inconsitency in the graph coloring, there is only one challenge edge that will reveal the cheating, and that inconsitency can easily be "shifted" around the graph to a random place at the commitment stage.

Let $E$ = number of edges

Let $R$ = number of rounds

The chance of the Prover getting away with cheating in a single round would be $\frac{E-1}{E}$

But every round is an independent event because the inconsitent edge is shifted, so after $R$ rounds the probability of not being caught is

$(\frac{E-1}{E})^R$

Let $E = 15, R = 15^2$, probability of cheating $= (\frac{14}{15})^{15^2} \approx 1.8 * 10^{-7}$

