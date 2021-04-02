import random
import math
import matplotlib.pyplot as plt
import numpy as np


# 编码
def geneEncoding(pop_size, chrom_length):
    pop = [[]]
    for i in range(pop_size):
        temp = []
        for j in range(chrom_length):
            temp.append(random.randint(0, 1))
        pop.append(temp)
    return pop[1:]  # ???


# 解码
def decodechrom(pop, chrom_length):
    temp = []
    for i in range(len(pop)):
        t = 0
        for j in range(chrom_length):
            t += pop[i][j] * (math.pow(2, chrom_length - 1 - j))
        temp.append(t)
    return temp


# 计算函数值
def calobjValue(pop, chrom_length, min_value, max_value):
    temp1 = []
    obj_value = []
    temp1 = decodechrom(pop, chrom_length)
    for i in range(len(temp1)):
        x = min_value + temp1[i] * max_value / (math.pow(2, chrom_length) - 1)
        obj_value.append(10 * math.sin(5 * x) + 7 * math.cos(4 * x))
    return obj_value


# 确定适应度函数
def calfitValue(obj_value, c_min):
    fit_value = []
    for i in range(len(obj_value)):
        if obj_value[i] < 0:
            temp = c_min - obj_value[i]
        else:
            temp = c_min
        fit_value.append(temp)
    return fit_value


# 找出最优解和最优解的基因编码
def best(pop, fit_value):
    px = len(pop)
    best_individual = []
    best_fit = fit_value[0]
    for i in range(1, px):
        if fit_value[i] > best_fit:
            best_fit = fit_value[i]
            best_individual = pop[i]
    return [best_individual, best_fit]


#
def sum(fit_value):
    total = 0
    for i in range(len(fit_value)):
        total += fit_value[i]
    return total


#
def cumsum(fit_value):
    for i in range(len(fit_value) - 2):
        t = 0
        j = 0
        while (j <= i):
            t += fit_value[j]
            j += 1
        fit_value[i] = t
        fit_value[len(fit_value) - 1] = 1


# 选择算子
def selection(pop, fit_value):
    newfit_value = []
    # 适应度总和
    total_fit = sum(fit_value)
    for i in range(len(fit_value)):
        newfit_value.append(fit_value[i] / total_fit)
    # 计算累积概率
    cumsum(newfit_value)
    ms = []
    pop_len = len(pop)
    for i in range(pop_len):
        ms.append(random.random())
    ms.sort()
    fitin = 0
    newin = 0
    newpop = pop
    # 轮盘赌选择法
    while newin < pop_len:
        if (ms[newin] < newfit_value[fitin]):
            newpop[newin] = pop[fitin]
            newin += 1
        else:
            fitin += 1
    pop = newpop


# 交配
def crossover(pop, pc):
    pop_len = len(pop)
    for i in range(pop_len - 1):
        if random.random() < pc:
            cpoint = random.randint(0, len(pop[0]))
            temp1 = []
            temp2 = []
            temp1.extend(pop[i][0:cpoint])
            temp1.extend(pop[i + 1][cpoint:len(pop[i])])
            temp2.extend(pop[i + 1][0:cpoint])
            temp2.extend(pop[i][cpoint:len(pop[i])])
            pop[i] = temp1
            pop[i + 1] = temp2


# 变异
def mutation(pop, pm):
    px = len(pop)
    py = len(pop[0])

    for i in range(px):
        if random.random() < pm:
            mpoint = random.randint(0, py - 1)
            if pop[i][mpoint] == 1:
                pop[i][mpoint] = 0
            else:
                pop[i][mpoint] = 1


# 计算2进制序列代表的（十进制）数值
def b2d(b, min_value, max_value, chrom_length):
    t = 0
    for j in range(len(b)):
        t += b[j] * (math.pow(2, len(b) - 1 - j))
    t = min_value + t * max_value / (math.pow(2, chrom_length) - 1)
    return t


# 找出最优解并在原函数图像上标出
def plot_optimal(results, c_min):
    best_fun = 0
    best_individual = 0
    for item in results:
        if item[0] > best_fun:
            best_fun = item[0]
            best_individual = item[1]
    best_fun = c_min - best_fun
    x = np.arange(0.0, 10.0, 0.001)
    y = []
    for item in x:
        y.append(10 * math.sin(5 * item) + 7 * math.cos(4 * item))
    plt.plot(x, y)
    plt.plot(best_individual, best_fun, 'r*')


if __name__ == '__main__':
    pop_size = 500  # 种群大小
    stop_generation = 1000  # 终止代数
    min_value = 0  # 定义域的左端点
    max_value = 10  # 定义域的右端点
    chrom_length = 10  # 染色体长度
    pc = 0.6  # 交叉概率
    pm = 0.001  # 变异概率
    results = [[]]  # 存储每一代的最优解，N个二元组
    fit_value = []  # 个体适应度
    fit_mean = []  # 平均适应度
    c_min = 1  # 适应度函数中的一个常数
    pop = geneEncoding(pop_size, chrom_length)

    for i in range(stop_generation):
        obj_value = calobjValue(pop, chrom_length, min_value, max_value)
        fit_value = calfitValue(obj_value, c_min)
        best_individual, best_fit = best(pop, fit_value)
        results.append([best_fit, b2d(best_individual, min_value, max_value,
                                      chrom_length)])
        selection(pop, fit_value)
        crossover(pop, pc)
        mutation(pop, pm)

    results = results[1:]
    # results.sort()

    X = []
    Y = []
    for i in range(stop_generation):
        X.append(i)
        Y.append(-results[i][0] + c_min)

    plt.plot(X, Y)
    plt.show()
    plot_optimal(results, c_min)
    plt.show()
