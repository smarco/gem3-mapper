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
mylogit <- glm(tp ~ edit+event+swg+sigma+mapq_1+mapq_2+sub_edit+sub_event+sub_swg+sub_sigma+ccand+acand+cmatch+mcs1+mcs2+arl+mrl+kmerf , data = data1, family = "binomial")

// SE
mylogit <- glm(tp ~ edit+event+swg+sub_edit+sub_event+sub_swg+ccand+cmatch+mcs1+mrl+kmerf , data = data1, family = "binomial")
// Uniq (chr1,HS)
// Mmaps (chr1,HS)
// Ties (chr1,HS)

/*
 * PE NEW
 */
// Uniq (chr1,HS)
// Mmaps (chr1,HS)
// Ties (chr1,HS)

