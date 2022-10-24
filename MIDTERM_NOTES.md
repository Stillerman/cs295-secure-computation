## Finite Fields
- It is not possible to pick a random integer because there are infinite
possibilities
- once you restrict yourself to a finite field, you can assign each number a
equal probability

## Galois Field
- Operations: +, -, *, division
- Why does the characteristic need to be prime?
    - division (multiplicative inverse) only exists if the characteristic is prime
- dont need p prime for additive secret sharing
- If p is not prime you have a Ring (not a Galois field)
- We need prime for Shamir secret sharing

## Ideal functionality
- Trusted third party who does all the computation
- Each party sends its stuff to the TTP
- TTP computes $F$
- Sends results back

## Formal Security
- If you can learn `x` in the real world, you can learn `x` in the ideal world

## Proving Security with Simulators Example

> Describe a *simulator* for the ideal world, that constructs views which are indistinguishable from the ones described in Question 1.

For party i we have
- Private input x_i
- final output of ideal func (sum of all xs)

We need to construct
- the input x_i (We have this explicitly from ideal func)
- shares of each other party's input (pick a random number because additive secret shares look like random numbers)
- shares of the final sum (we have this as output from ideal func)

Because we cannot distinguish between real random numbers and the additive secret random numbers, we cannot distinguish from real life vs protocol and we are secure.

## Types of secret sharing
|  | Additive | Shamir |
|--|----|---|
|format| 1 field elm. | (x,y) coord of field elm. |
|reconstruction|sum|polyn. interpolation|
|homomorphism|additive|additive and multiplicative (w/ degree reduction)|
|shares needed to reconstruct|all n|t (threshold)|
|use cases|Binary GF(2)|GF(p) for large prime p.|

## Calculating Mean (additive)

**Protocol: Secure Summation with Additive Secret Sharing**
- **Round 1**: Each party $P_i$ sends one share of its input $x_i$ to each other party, keeping one share for itself.
- **Round 2**: Each party $P_i$ sums the shares it holds (including both the shares it has received and the share it kept for itself). Each party sends its sum to all other parties.
- **Opening**: Each party adds up the sums it receives and the sum it computed and then divides by the quantity of sums recieved. The quantity of sums recieves represents the number of parties which is equal to $n$.

```py
class MeanParty(Party):
    def round1(self, parties, input_num):
        self.input = GF(input_num)
        self.parties = parties
        n = len(parties)
        t = n-1
        
        shares = shamir_share(self.input, t, n)

        for share, party in zip(shares, parties):
            self.send(party, 1, share)

    def round2(self):
        baseline_x_coord = self.received[1][0][0]
        cummulative_y_coord = GF(0)

        for xcoord, ycoord in self.received[1]:
            assert xcoord == baseline_x_coord
            cummulative_y_coord += ycoord

        for party in self.parties:
            self.send(party, 2, (baseline_x_coord, cummulative_y_coord))

    
    def round3(self):
        # reconstruct the sum
        sum_share = reconstruct(self.received[2])

        self.output = int(sum_share) / len(self.parties)
```

## Multiply three numbers (degree reduction)

- R1
  - Each party $P_i$ recieves shares $a_i$, $b_i$, $c_i$ as input
  - Let $s_i = a_i * b_i$ - threshold for $s_i$ will be less than or equal to $2t$ where t is the initial threshold
    - $s_i$ is exactly $q(\alpha_i)$
  - $P_i$ computes $h_i^1 ... h_i^n$ = `share`($s_i$, t, n)
  - $P_i$ sends share $h_i^j$ to party $j$
- R2
  - Each party $P_i$ recieves the shares $h_j^i$ (yes the sub and superscripts are flipped)
  - $P_i$ computes $\sum_j (h_j^i * \lambda_j)$. This value is the product of $a$ and $b$ but not yet $c$. So then this product is shamir shared with threshold $t$, and similarly to round 1, shares of $c$ are multiplied in and then degree reduction started.
- R3
  - Each party $P_i$ recieves the shares $h_j^i$ (yes the sub and superscripts are flipped)
  - $P_i$ computes $\sum_j (h_j^i * \lambda_j)$. This is the final product

```py
class MultThreeParty(Party):
    def round1(self, parties, a_shr, b_shr, c_shr, t):
        self.input = (a_shr, b_shr, c_shr)
        self.c_shr = c_shr # save this one for later
        self.parties = parties
        n = len(parties)
        assert t <= n/2
        
        # - Each party $P_i$ recieves shares $a_i$, $b_i$ as input
        # - Let $s_i = a_i * b_i$ - threshold for $s_i$ will be less than or equal to $2t$ where t is the initial threshold
        #     - $s_i$ is exactly $q(\alpha_i)$

        a_x, a_y = a_shr
        b_x, b_y = b_shr
        c_x, c_y = c_shr

        # they better have the same x coord
        assert a_x == b_x == c_x

        s_i = a_y * b_y #q(x_i) higher degree than we'd like (degree 2t at most)
        self.x_coord = a_x # save this for round 2

        # - $P_i$ computes $h_i^1 ... h_i^n$ = `share`($s_i$, t, n)
        h_i_js = shamir_share(s_i, t, n)
        # - $P_i$ sends share $h_i^j$ to party $j$
        for party, share in zip(self.parties, h_i_js):
            self.send(party, 1, share)

    def round2(self):
        n = len(self.parties)
        
        # - Each party $P_i$ recieves the shares $h_j^i$ (yes the sub and superscripts are flipped)
        h_j_is = self.received[1]
        h_j_is_y = [s[1] for s in h_j_is]

        # $P_i$ computes $\sum_j (h_j^i * \lambda_j)$ and ouptuts this value as its own share of the origional product with threshold $t$
        V_a = GF(np.vander(range(1,n+1), increasing=True))
        V_a_inv = np.linalg.inv(V_a)
        lambda_js = V_a_inv[0]

        prods = [h_j_is_y[i] * lambda_js[i] for i in range(n)]

        self.new_share = self.x_coord, GF(prods).sum()

        n_x, n_y = self.new_share
        c_x, c_y = self.c_shr

        s_i = n_y * c_y

        h_i_js = shamir_share(s_i, t, n)

        for party, share in zip(self.parties, h_i_js):
            self.send(party, 2, share)

    def round3(self):
        n = len(self.parties)

        h_j_is = self.received[2]
        h_j_is_y = [s[1] for s in h_j_is]

        # $P_i$ computes $\sum_j (h_j^i * \lambda_j)$ and ouptuts this value as its own share of the origional product with threshold $t$
        V_a = GF(np.vander(range(1,n+1), increasing=True))
        V_a_inv = np.linalg.inv(V_a)
        lambda_js = V_a_inv[0]

        prods = [h_j_is_y[i] * lambda_js[i] for i in range(n)]

        self.output = self.x_coord, GF(prods).sum()

        return self.output
```
### Test Case
```py
# TEST CASE for question 2

NUM_PARTIES = 6
# (t, n)-Shamir scheme
n = NUM_PARTIES
t = 3

shares1 = shamir_share(5, t, n)
shares2 = shamir_share(6, t, n)
shares3 = shamir_share(7, t, n)

parties = [MultThreeParty() for _ in range(NUM_PARTIES)]

for p,s1,s2,s3 in zip(parties, shares1, shares2, shares3):
    p.round1(parties, s1, s2, s3, t)
for p in parties:
    p.round2()
for p in parties:
    p.round3()
for p in parties:
    # print(p.get_view())
    print(p.output)

output_shares = [p.output for p in parties]
print('Reconstruction, with all shares:', reconstruct(output_shares))
print('Reconstruction, with 3 shares:', reconstruct(output_shares[:3]))
print('Reconstruction, with 2 shares:', reconstruct(output_shares[:2]))

assert reconstruct(output_shares) == 210
assert reconstruct(output_shares[:3]) == 210
assert reconstruct(output_shares[:2]) != 210
```

## BGW

- Round 1 each party P_i generates shamir share of each of its secret inputs, sends one share to each party
- Round 2 each party P_i receives shares of its inputs from each party and initializes the wire_vals dict
- Round n
    - evaluarte the next gate in the circuit
    - if add, add the values of the two inputs and store in the output wire
    - if mult, multiply the values of the two inputs and store in the output wire, perform degree reduction
- #round n+1
    - each party P_i broadcasts its shares of output wires
- round n+2
    - each party P_i reconstructs the output wires from the shares it received

```py
class BGWParty(Party):
    def round1(self, parties, circuit, my_inputs):
        self.parties = parties
        self.is_done = False
        self.circuit = circuit
        n = len(parties)
        t = int(n/2)

        # Round 1 (phase 1): Each party P_i create n Shamir shares of ech of its secret inputs ,
        # and sends one share to each other party
        my_id = parties.index(self)
        my_input_wires = circuit.inputs[my_id]
        # print(f"party num {my_id} my input wires {my_inputs}")

        # input shares will map each party to the shares of my inputs destined for that party
        input_shares = {p: {} for p in parties}

        for wire, value in zip(my_input_wires, my_inputs):
            shares = shamir_share(value, t, n)
            for p, s in zip(parties, shares):
                input_shares[p][wire] = s

        for p in parties:
            self.send(p, 1, input_shares[p])


    def round2(self, my_id):
        self.wire_vals = {}

        # Round 2 (phase 2): EAch party recieves one share for
        # each input wire and initializs the wire_vals dict
        received_shares = self.received[1]

        for received_dict in received_shares:
            for key, value in received_dict.items():
                self.wire_vals[key] = value

        self.phase = 3
        self.current_gate = 0
        self.need_degree_reduction = False

    def roundn(self, round_num):
        n = len(self.parties)
        t = int(n/2)

        if self.need_degree_reduction:
            # finish the degree reduction and
            # update wire_vals
            # - each party $P_i$ receives shares $h_j^i$
            h_j_is = self.received[round_num - 1]
            h_js_is_y = [s[1] for s in h_j_is]

            V_a = GF(np.vander(range(1, n+1), increasing=True))
            V_a_inv = np.linalg.inv(V_a)
            lambda_js = V_a_inv[0] # first row

            prods = [lambda_j * s for lambda_j, s in zip(lambda_js, h_js_is_y)]

            #
            g = self.circuit.gates[self.current_gate]
            self.wire_vals[g.out] = \
                (self.x_coord, GF(prods).sum())

            self.current_gate += 1
            self.need_degree_reduction = False

        if self.current_gate >= len(self.circuit.gates) and self.phase == 3:
            self.phase = 4

        if self.phase == 3:
            # Evaluate the next gate in the circuit
            # If it is an ADD gate, look up the shares of its input sin the dict and add them together,then update the dict to map its output to the resulting share
            # If it is a MULT gate then, look up the shares of its input sin the dict and multiply them together,then perform degree reduction
            g = self.circuit.gates[self.current_gate]


            x1, y1 = self.wire_vals[g.in1] # lookup the value of the first input
            x2, y2 = self.wire_vals[g.in2]
            assert x1 == x2

            if g.type == 'ADD':
                self.wire_vals[g.out] = (x1, y1 + y2)
                self.current_gate += 1

            elif g.type == 'MULT':
                mult_result = y1 * y2
                # remember this value
                # setup the degree reduction
                # next time I enter `roundn` function,
                # finish the degree reduction and
                # update wire_vals
                self.x_coord = x1
                self.need_degree_reduction = True

                h_i_js = shamir_share(mult_result, t, n)
                for party, share, in zip(self.parties, h_i_js):
                    self.send(party, round_num, share)


        elif self.phase == 4:
            # Round k: (phase 4) When all gates have been evaluated
            # Each party P_u broadcasts its shares of output wires
            output_wires = self.circuit.outputs
            output_shares = [self.wire_vals[w] for w in output_wires]

            for p in self.parties:
                self.send(p, round_num, output_shares)

            self.phase = 5

        elif self.phase == 5:
            # Round k+1: (phase 5) each party receives n shares of each output wire value, reconstructs each wire’s actual value, and outputs the values
            received_shares = self.received[round_num - 1]

            output_shares = [ [] for _ in self.circuit.outputs]
            # arrange the shares
            for shares in received_shares:
                # shares received from a single party p_i
                for j, wire_share in enumerate(shares):
                    # this is the share for wire j
                    output_shares[j].append(wire_share)

            # do the reconstruction
            output_vals = []
            for shares in output_shares:
                output_vals.append(reconstruct(shares))

            self.output = output_vals

            self.is_done = True
```

## OT
Describe a method for evaluating an `AND` gate using 1-out-of-4 OT on additive-secret-shared inputs.

- P1 and P2 each hold one additive share of the two input wire values
- We use 1 out of 4 OT with P1 as S and P2 as R
- P1 will generate their output share randomly
- We will build a truth table for P2's output share
- P1 can say
  - Given my input shares and my random output shares
  - What would p2's output share be for each possible value of P2's input shares
  - We build a table of these, and use the potential output shares for P2 as the OT secrets
  - P2 uses its shares as the OT selection bits


**Protocol**:
- Inputs: P1 has s1_i, s1_j; P2 has s2_i, s2_j (additive shares of the AND gates input)


**Round1** P2 generates keypaires and keep one of them
- P2 generates 4 keypairs
- P2 keeps one secret key based on the values of s2_i and s2_j
  - if s2_i = 0 and s2_j = 0, R sends (pk1, _, _, _) to S
  - ...


**Round2** P1 generates the truth table as its secrets and encrypts the values of it
- P1 recieves 4 public keys from P2
- P1 generates a random output share r = s1_k
- P1 calls T_G to get the truth table, using s1_i, s1_j, and r as inputs
- P1 encrypts each row of the truth table using the 4 public keys


**Round3** P2 decrypts the row of the truth table corresponding to its actual shares
- P2 decrypts the right row of the trith table using sk1
  
At the end, P1 has s1_k, P2 has s2_k and s1_k + s2_k = output of the gate

```py
class AND_P1(Party):
    # x1 and x2 are the secrets
    def round1(self, s1_i, s1_j, p2):
        self.s1_i = s1_i
        self.s1_j = s1_j
        self.p2 = p2

    def round2(self):
        # Round 2: S receive (pka, pkb). S sends (Enc_pka(x1), Enc_pkb(x2)) to R
        [pks] = self.received[1]
       
#         P1 generates a random output share r = s1_k
#         P1 calls T_G to get the truth table, using s1_i, s1_j, and r as inputs
        r = GF_2.Random()
        self.output = r
        truth_table = T_G(r, self.s1_i, self.s1_j)
        encrypted_truth_table = []
        for pk, table_element in zip(pks, truth_table):
            table_element_b = int(table_element).to_bytes(1, 'little')
            enc = SealedBox(pk).encrypt(table_element_b)
            encrypted_truth_table.append(enc)
       
        self.send(self.p2, 2, encrypted_truth_table)
   
    def round3(self):
        return self.output

class AND_P2(Party):
    def round1(self, s2_i, s2_j, p1):
        self.p1 = p1
        self.s2_i = s2_i
        self.s2_j = s2_j
#         P2 generates 4 keypairs
#         P2 keeps one secret key, based on the values of s2_i, s2_j
        # Round 1: R generates two keypairs: sk1, pk1 and sk2, pk2. R throws away sk2.
        keypair1 = PrivateKey.generate() # keep this one
        keypair2 = PrivateKey.generate() # throw this one away after this round
        keypair3 = PrivateKey.generate() # throw this one away after this round
        keypair4 = PrivateKey.generate() # throw this one away after this round

        self.saved_key = keypair1
       
        if s2_i == 0 and s2_j == 0:
            self.send(self.p1, 1, (keypair1.public_key,
                                   keypair2.public_key,
                                   keypair3.public_key,
                                   keypair4.public_key))
        elif s2_i == 0 and s2_j == 1:
            self.send(self.p1, 1, (keypair2.public_key,
                                   keypair1.public_key,
                                   keypair3.public_key,
                                   keypair4.public_key))
        elif s2_i == 1 and s2_j == 0:
            self.send(self.p1, 1, (keypair3.public_key,
                                   keypair2.public_key,
                                   keypair1.public_key,
                                   keypair4.public_key))
        elif s2_i == 1 and s2_j == 1:
            self.send(self.p1, 1, (keypair4.public_key,
                                   keypair2.public_key,
                                   keypair3.public_key,
                                   keypair1.public_key))

   
    def round2(self):
        pass
   
    def round3(self):
        [(c1, c2, c3, c4)] = self.received[2]
       
        if self.s2_i == 0 and self.s2_j == 0:
            plaintext = SealedBox(self.saved_key).decrypt(c1)
        elif self.s2_i == 0 and self.s2_j == 1:
            plaintext = SealedBox(self.saved_key).decrypt(c2)
        elif self.s2_i == 1 and self.s2_j == 0:
            plaintext = SealedBox(self.saved_key).decrypt(c3)
        elif self.s2_i == 1 and self.s2_j == 1:
            plaintext = SealedBox(self.saved_key).decrypt(c4)
       
        self.output = GF_2(int.from_bytes(plaintext, 'little'))
        return self.output
```

## GMW

* Refer to 10-03-22 Q3 to write GMW *

**GMW**
- Additive Shares
- Eval Gates
  - AND gate uses OT
- Reconstruct outputs


**Inputs**
- P1 has (binary) values for some of the input wires
- P2 has (binary) values for some of the input wires

**Protocol**
- Input section
  - **Round 1** (phase 1): Each party secret shares its input bits to the other party (same as BGW but with additive) and adds its own shares to the wire_vals dictionary (create two shares and put one in dict and send other to other party)
  - **Round 2** each party recieves shares of the other parties inputs and adds these to the wire_vals dictionary


- Eval section
  - ** R 2 < i < k** Eval next gate $g$ in circuit
    - XOR: `wire_vals[g.out] = wire_vals[g.in1] + wire_vals[g.in2]` (similar to addition gate in BGW)
    - INV (NOT): `wire_vals[g.out] = wire_vals[g.in1] + GF_2(1)`
      - `g.in2 = -1` in this case to indicate that its not used
    - AND: hard case - use the protocol from Q7
      - inputs: wire_vals[g.in1] and wire_val[g.in2]
      - outputs: one secret share for g.out
      - implement 3 ot-phases
        - OTP1: generate pub keys (r1 from Q7)
        - OTP2: generate and encrypt truth table (r2)
        - OPT3: decrypt one row from the TT (r3)
- Output section
  - ** round k ** parties broadcast their shares of the output wire values
  - ** round k+1** parties reconstruct the output wire values using their own shares plus the broadcasted shares recieved from the other party

## Circuit based protocol comparison
Circuit has gates
- Binary (AND, OR, NOT)
- Arithmetic circuits (*, +)
- Strengths and weaknesses of each
  - Binary Strengths:
    - Boolean functions
    - Equality
    - Comparisons(Arithmetic bad at this)
  - Arithmetic Strengths:
    - Smaller
    - Faster

| |BGW|GMW|
|--|--|--|
|circuits|arithmetic|binary|
|format|shamir shares|additive|
|parties|$n$|2|
|round complexity|$O(d)$|$O(d)$|