from collections import deque, namedtuple

Minimizer = namedtuple("Minimizer", ["kmer", "pos"])


def minimize_func(k, w, hash_func):
    def minimize(sequence):
        hashes = deque([hash_func(sequence[i:i+k]) for i in range(w-k+1)])
        cur_min = min(hashes)
        min_idxs = [i for i, h in enumerate(hashes) if h == cur_min]
        minimizers = [Minimizer(cur_min, i) for i in min_idxs]
        cur_min_count = len(min_idxs)
        for i in range(w-k+1, len(sequence)-k+1):
            print("#", hashes)
            cur_hash = hashes.popleft()
            if cur_hash == cur_min:
                cur_min_count -= 1
            new_hash = hash_func(sequence[i:i+k])
            hashes.append(new_hash)
            # print(i, new_hash, cur_hash)
            if new_hash < cur_min:
                cur_min = new_hash
                cur_min_count = 1
                minimizers.append(Minimizer(cur_min, i))
            elif new_hash == cur_min:
                cur_min_count += 1
                minimizers.append(Minimizer(cur_min, i))
            elif cur_min_count == 0:
                cur_min = min(hashes)
                min_idxs = [j+i-(w-k) for j, h in enumerate(hashes) if h == cur_min]
                minimizers.extend(Minimizer(cur_min, j) for j in min_idxs)
                cur_min_count += 1
        minimizers.sort()
        return minimizers

    return minimize


if __name__ == "__main__":
    f = minimize_func(2, 4, lambda x: x)
    print(f([1, 2, 3, 0, 1, 2, 2, 3, 1]))
