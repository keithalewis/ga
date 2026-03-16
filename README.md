# Grassmann Algebra

Hermann Grassmann invented an algebra for Euclidean geometry.

Let $E$ be points in space. His only rule for (exterior) products is that
they are associative and $PQ = 0$ if and only if $P = Q$.
We have $PQ = -PQ$ since $0 = (P + Q)(P + Q) = PQ + QP$.

The _order_ of $E$ is the largest $n$ for which there exist $P_j\in E$
with $P_0 P_1 \cdots P_n \not= 0$. Whenever this holds we say the points
are _independent_. The _dimension_ is one less than the order.

A fundamental results is that if $P\in E$ and $P P_0 \cdots P_k = 0$ and
$P_0,\ldots P_k$ are independent then

$P = \sum_{j = 0}^k P_0 \cots P_{j-1} P P_{j + 1}\cdots P_k/P_0\cdots P_k$.

$P$ belongs to the smallest convex set containing $P_0,\ldots, P_k$ if
and only if all coefficients are non-negative. $P$ belongs to the smallest
cone with vertex $P_0$ is and only if all coefficients except the
first are non-negative.