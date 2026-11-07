https://github.com/vvvooolt/Attaching-edges-in-the-supersingular-isogeny-spine# Verification code for "Attaching edges in the supersingular isogeny spine"

Four SageMath files check the paper's results against isogeny graphs built from scratch. The paper's
theorems are proved; these checks are an independent test of the statements and of the computed facts.

| File              | Checks                                                                                                                                                                                                                                       |
| ----------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `lib.sage`        | shared helpers: curves and isogenies over F_p (ground truth), binary quadratic forms (predictions), class-group action on the floor                                                                                                          |
| `identities.sage` | every fact marked "computed" in the paper: H₋₁₅ (Lemma 3.1), the resultants of Table 1 for ℓ = 2, 3, 5, 7, Remark 3.12, Table 2, the ℓ = 3 data of Section 6.2, the H₋₁₄₀ identity (Cor 6.3), and the numbers in Example 1.1 and Section 3.8 |
| `check_l2.sage`   | ℓ = 2: the edge is new, Lemma 3.7, floor (Cor 3.6), Theorem A, Theorem C with the exact length (Thm 3.17) and Cor 3.18, and Theorem B by labelling every floor curve with its class                                                          |
| `check_odd.sage`  | odd ℓ: Cor 6.1 (the five prime-degree families of Table 2), Cor 6.4 (D′ = −20), Cor 6.6 (D′ = −32): levels, attachment, exact lengths, and the bounds of Theorem E                                                                           |

**Ground truth** uses only curves over F_p: the spine component is found by breadth-first search over the
F_p-rational ℓ-isogenies of both twists of each j-invariant, and "new" means a double root of Φ_ℓ(j₁, Y)
with no rational ℓ-isogeny behind it. **Predictions** use only binary quadratic forms. The two never share
code.

## Running

    sage identities.sage                 # under a minute
    sage check_l2.sage 60 100000 5000    # Theorem B checked for p < 5000; a few hours on one core
    sage check_odd.sage 0 50000          # optional: family names, e.g. 3,-35 7,-187 (hours on one core)

Each prints one line per prime and ends with a `SUMMARY` line listing any failures. Long ranges can be split
into chunks and run in parallel.

## Results (`results/`)

| Run              | Range                                                            | Primes                                                | Failures |
| ---------------- | ---------------------------------------------------------------- | ----------------------------------------------------- | -------- |
| `identities.txt` | —                                                                | 26 checks                                             | 0        |
| `l2.txt`         | all p ≡ 11, 14 mod 15, 60 < p < 100,000; Theorem B for p < 5,000 | 2,392 (1,877 attaching, 515 not)                      | 0        |
| `odd.txt`        | the seven families, above each family's bound, p < 50,000        | 8,381 family–prime pairs (5,032 attaching, 3,349 not) | 0        |

## Sealed predictions (`sealed/`)

The paper reports four batches whose predictions were sealed before the ground truth was computed. Each
manifest lists the SHA-256 hashes of the prediction files and of the scripts that produced them (those
original scripts are included here only for that reason). The `.ots` file is an OpenTimestamps proof that
the manifest existed at the time of the Bitcoin block below; `*_compare.txt` is the comparison run after
confirmation.

| Paper                                                                      | Manifest                      | Bitcoin block |
| -------------------------------------------------------------------------- | ----------------------------- | ------------- |
| Cor 6.2 (D′ = −35)                                                         | `logs/v5b_manifest.txt`       | 969151        |
| Cor 6.6 (D′ = −32), sealed before the proof                                | `logs/v5d_manifest.txt`       | 969156        |
| Cor 6.1 (prime degrees)                                                    | `logs/v5e_manifest.txt`       | 969164        |
| Theorem D, ℓ ≤ 13                                                          | `logs/v7_oos/manifest.txt`    | 969197        |
| Cor 6.4 (D′ = −20): time stamp later than the run, so _not_ pre-registered | `logs/v5c_manifest.txt`       | 969156        |
| ℓ = 2 out-of-sample: hash recorded first, no blockchain stamp              | `logs/oos_predictions.sha256` | —             |

Check all hashes, from inside `sealed/`:

    for m in logs/v5b_manifest.txt logs/v5c_manifest.txt logs/v5d_manifest.txt logs/v5e_manifest.txt logs/v7_oos/manifest.txt; do
      sha256sum -c <(awk 'NF==2' "$m"); done
    sha256sum -c <(head -1 logs/oos_predictions.sha256)

Verifying the `.ots` proofs themselves needs the OpenTimestamps client (`ots verify`) and access to a Bitcoin node.
