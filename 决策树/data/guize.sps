STRING pre_001 (A4).

/* Node 7 */.
DO IF (是否周末 EQ "是")  AND  (天气 NE "好")  AND  (是否有促销 NE "否").
COMPUTE nod_001 = 7.
COMPUTE pre_001 = '高'.
COMPUTE prb_001 = 0.800000.
END IF.
EXECUTE.

/* Node 8 */.
DO IF (是否周末 EQ "是")  AND  (天气 NE "好")  AND  (是否有促销 EQ "否").
COMPUTE nod_001 = 8.
COMPUTE pre_001 = '低'.
COMPUTE prb_001 = 0.666667.
END IF.
EXECUTE.

/* Node 4 */.
DO IF (是否周末 EQ "是")  AND  (天气 EQ "好").
COMPUTE nod_001 = 4.
COMPUTE pre_001 = '高'.
COMPUTE prb_001 = 1.000000.
END IF.
EXECUTE.

/* Node 9 */.
DO IF (是否周末 NE "是")  AND  (是否有促销 NE "否")  AND  (天气 EQ "坏").
COMPUTE nod_001 = 9.
COMPUTE pre_001 = '低'.
COMPUTE prb_001 = 0.600000.
END IF.
EXECUTE.

/* Node 10 */.
DO IF (是否周末 NE "是")  AND  (是否有促销 NE "否")  AND  (天气 NE "坏").
COMPUTE nod_001 = 10.
COMPUTE pre_001 = '高'.
COMPUTE prb_001 = 0.571429.
END IF.
EXECUTE.

/* Node 11 */.
DO IF (是否周末 NE "是")  AND  (是否有促销 EQ "否")  AND  (天气 EQ "坏").
COMPUTE nod_001 = 11.
COMPUTE pre_001 = '低'.
COMPUTE prb_001 = 1.000000.
END IF.
EXECUTE.

/* Node 12 */.
DO IF (是否周末 NE "是")  AND  (是否有促销 EQ "否")  AND  (天气 NE "坏").
COMPUTE nod_001 = 12.
COMPUTE pre_001 = '低'.
COMPUTE prb_001 = 0.750000.
END IF.
EXECUTE.
