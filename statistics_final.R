# 初步工作
library(ISLR);
data=Auto[,-9];
dataYear=split(data,data$year);
dataOrigin=split(data,data$origin);
dataCylinders=split(data,data$cylinders);

# 均值
print(mean(data$mpg));
print(mean(data$displacement));
print(mean(data$horsepower));
print(mean(data$weight));
print(mean(data$acceleration));

# 中位数
print(median(data$mpg));
print(median(data$displacement));
print(median(data$horsepower));
print(median(data$weight));
print(median(data$acceleration));

# 极差
print(max(data$mpg)-min(data$mpg));
print(max(data$displacement)-min(data$displacement));
print(max(data$horsepower)-min(data$horsepower));
print(max(data$weight)-min(data$weight));
print(max(data$acceleration)-min(data$acceleration));

# 修正样本方差
print(var(data$mpg));
print(var(data$displacement));
print(var(data$horsepower));
print(var(data$weight));
print(var(data$acceleration));

# 箱线图
library(ggplot2);
ggplot(data,aes(y=mpg))+geom_boxplot()
ggplot(data,aes(y=displacement))+geom_boxplot()
ggplot(data,aes(y=horsepower))+geom_boxplot()
ggplot(data,aes(y=weight))+geom_boxplot()
ggplot(data,aes(y=acceleration))+geom_boxplot()

# Q-Q图
qqnorm(data$mpg)
qqline(data$mpg)
qqnorm(data$acceleration)
qqline(data$acceleration)

# 偏度峰度
library(moments)
skewness(data$mpg)
kurtosis((data$mpg))
skewness(data$acceleration)
kurtosis((data$acceleration))
#skewness(data$horsepower)
#kurtosis((data$horsepower))

# 列出不同年、不同气缸数车型数量
table(dataCylinders$"4"$year)
table(dataCylinders$"4"$origin)
table(dataCylinders$"6"$year)
table(dataCylinders$"6"$origin)
table(dataCylinders$"8"$year)
table(dataCylinders$"8"$origin)

# 独立性检验：手算

# 斯米尔诺夫检验
table(data$year)
ks.test(dataYear$"72"$horsepower,dataYear$"77"$horsepower)
# print(sort(dataYear$"72"$horsepower))
# print(sort(dataYear$"77"$horsepower))
# wilcox.test(dataYear$"72"$horsepower,dataYear$"77"$horsepower)

# 单因素方差分析
n4Cyl=nrow(dataCylinders$"4");
n6Cyl=nrow(dataCylinders$"6");
n8Cyl=nrow(dataCylinders$"8");
nCyl=n4Cyl+n6Cyl+n8Cyl;
var4Cyl=var(dataCylinders$"4"$mpg)*(n4Cyl-1)/n4Cyl;
var6Cyl=var(dataCylinders$"6"$mpg)*(n6Cyl-1)/n6Cyl;
var8Cyl=var(dataCylinders$"8"$mpg)*(n8Cyl-1)/n8Cyl;
ita4Bar=mean(dataCylinders$"4"$mpg);
ita6Bar=mean(dataCylinders$"6"$mpg);
ita8Bar=mean(dataCylinders$"8"$mpg);
itaBar=(ita4Bar*n4Cyl+ita6Bar*n6Cyl+ita8Bar*n8Cyl)/nCyl;
QA=ita4Bar*ita4Bar*n4Cyl+ita6Bar*ita6Bar*n6Cyl;
QA=QA+ita8Bar*ita8Bar*n8Cyl-itaBar*itaBar*nCyl;
QE=n4Cyl*var4Cyl+n6Cyl*var6Cyl+n8Cyl*var8Cyl;
print((nCyl-2)*QA/2/QE);

ita3Bar=mean(dataCylinders$"3"$mpg);
ita5Bar=mean(dataCylinders$"5"$mpg);
# 线性回归分析：只有四个样本，手算即可


