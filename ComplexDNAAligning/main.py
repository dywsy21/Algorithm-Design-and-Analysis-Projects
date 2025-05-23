from evaluate import calculate_value

def find_alignment_tuples(query, ref):
    pass


def main():
    q1 = q2 = r1 = r2 = ''
    with open('data/query1.txt', 'r') as f:
        q1 = f.read()
    with open('data/query2.txt', 'r') as f:
        q2 = f.read()
    with open('data/ref1.txt', 'r') as f:
        r1 = f.read()
    with open('data/ref2.txt', 'r') as f:
        r2 = f.read()
    q1 = q1.strip()
    q2 = q2.strip()
    r1 = r1.strip()
    r2 = r2.strip()
    t1 = find_alignment_tuples(q1, r1)
    t2 = find_alignment_tuples(q2, r2)
    print('t1: ', t1)
    print('t2: ', t2)
    print('t1 score: ', calculate_value(t1, r1, q1))
    print('t2 score: ', calculate_value(t2, r2, q2))
