rm(list=ls())
data1=read.table('http://quantedu.com/wp-content/uploads/2014/04/SAT.txt',header=T)
head(data1,2)
data1$Race=factor(data1$Race)
mylogit <- glm(Admit ~ SAT + GPA + Race, data = data1, family = "binomial")
summary(mylogit)
newdata = data.frame(SAT=1600, GPA=4, Race=factor(1))
predict(mylogit, newdata, type="response") 
newdata1 = data.frame(SAT=1600, GPA=4, Race=factor(2))
predict(mylogit, newdata1, type="response") 
exp(cbind(OR = coef(mylogit), confint(mylogit)))

// Uniq
rm(list=ls())
data1=read.table('logit.uniq',header=T)
mylogit <- glm(tp ~ event+swg+mcs, data = data1, family = "binomial")
mylogit
summary(mylogit)
exp(cbind(OR = coef(mylogit), confint(mylogit)))

/*
 * SE NEW
 */
// ALL
mylogit <- glm(tp ~ edit+event+swg+sigma+mapq_1+mapq_2+edit1+edit2+event1+event2+swg1+swg2+sigma1+sigma2+ack+mcs1+mcs2+mrl+kmerf , data = data1, family = "binomial")

// SE
mylogit <- glm(tp ~ edit+event+swg+edit1+edit2+event1+event2+swg1+swg2+ack+mcs1+mrl+kmerf , data = data1, family = "binomial")
// Uniq (chr1,HS)
mylogit <- glm(tp ~ edit+event+mcs1+mrl+kmerf , data = data1, family = "binomial")
// Mmaps (chr1,HS)
mylogit <- glm(tp ~ edit+event+event2+swg1+swg2+mcs1+mrl , data = data1, family = "binomial")
// Ties (chr1,HS)
mylogit <- glm(tp ~ edit+event+event2+mcs1+mrl, data = data1, family = "binomial")

/*
 * PE NEW
 */
// Uniq (chr1,HS)
mylogit <- glm(tp ~ edit+event+mcs1+mcs2+kmerf, data = data1, family = "binomial")
mylogit <- glm(tp ~ edit+event+mcs1+mcs2+kmerf , data = data1, family = "binomial")
// Mmaps (chr1,HS)
mylogit <- glm(tp ~ edit+event+swg+sigma+mapq_1+mapq_2+edit1+edit2+event2+swg1+swg2+sigma1+sigma2+ack+mcs1+mcs2+mrl+kmerf , data = data1, family = "binomial")
// Ties (chr1,HS)
mylogit <- glm(tp ~ edit+event+swg+sigma+mapq_1+mapq_2+sigma1+ack+mcs1+mcs2+mrl, data = data1, family = "binomial")
