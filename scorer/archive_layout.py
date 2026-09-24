"""Strict layout adapter for archived variable-atom-count candidate batches."""
import math
from numbers import Integral


def shard_indices(total, batch_size, num_shards, shard_idx):
    if any(not isinstance(x, Integral) for x in (total, batch_size, num_shards, shard_idx)):
        raise ValueError('Integer shard parameters required')
    if total < 0 or batch_size < 1 or num_shards < 1 or not 0 <= shard_idx < num_shards:
        raise ValueError('Invalid shard parameters')
    return [i for start in range(0, total, batch_size)
            if (start // batch_size) % num_shards == shard_idx
            for i in range(start, min(start + batch_size, total))]


def unpack_spans(counts, batch_size, packed_length):
    if not counts or not counts[0] or batch_size < 1 or packed_length < 0:
        raise ValueError('Nonempty candidate-count matrix and valid dimensions required')
    n = len(counts[0])
    if any(len(row) != n for row in counts):
        raise ValueError('Ragged candidate-count matrix')
    if any(not isinstance(v, Integral) or v < 0 for row in counts for v in row):
        raise ValueError('Atom counts must be nonnegative integers')
    spans = [[] for _ in counts]
    offset = 0
    for start in range(0, n, batch_size):
        end = min(start + batch_size, n)
        width = max(sum(row[start:end]) for row in counts)
        for k, row in enumerate(counts):
            cursor = offset
            for value in row[start:end]:
                spans[k].append((cursor, cursor + value))
                cursor += value
        offset += width
    if offset != packed_length:
        raise ValueError(f'Packed length {packed_length} differs from required {offset}')
    return spans


def summarize_candidates(distances, budgets=(1, 5, 10, 20)):
    if any(not isinstance(k, Integral) or k < 1 or k > len(distances) for k in budgets):
        raise ValueError('Incomplete or invalid candidate budget')
    if any(d is not None and (not math.isfinite(d) or d < 0) for d in distances):
        raise ValueError('Invalid RMSD value')
    result = []
    for k in budgets:
        hits = [d for d in distances[:k] if d is not None]
        result.append({'k': k, 'matched': bool(hits), 'best_rmsd': min(hits) if hits else None})
    return result
