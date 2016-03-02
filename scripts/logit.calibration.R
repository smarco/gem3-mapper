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
// None
mylogit <- glm(tp ~ , data = data1, family = "binomial")
// ALL
mylogit <- glm(tp ~ edit+sub_edit+event+sub_event+swg+sub_swg+mcs_1+fs_matches+ack_1+max_region+mp+m2p, data = data1, family = "binomial")

// Uniq (chr1,HS)
mylogit <- glm(tp ~ edit+event+swg+mcs_1+max_region+mp, data = data1, family = "binomial")
// Mmaps (chr1,HS)
mylogit <- glm(tp ~ edit+event+swg+sub_swg+mcs_1+max_region+mp, data = data1, family = "binomial")
// Ties (chr1,HS)
mylogit <- glm(tp ~ edit+sub_edit+event+mcs_1+mp, data = data1, family = "binomial")

/*
 * PE NEW
 */
rm(list=ls())
data1=read.table('logit.uniq',header=T)
summary(mylogit)
exp(cbind(OR = coef(mylogit), confint(mylogit)))

// None
mylogit <- glm(tp ~ , data = data1, family = "binomial")
// ALL
mylogit <- glm(tp ~ edit+sub_edit+event+sub_event+swg+sub_swg+fs_matches+mcs_1+mcs_2+ack_1+ack_2+max_region+mp+m2p+sigma+sub_sigma+mapq_1+mapq_2, data = data1, family = "binomial")

// Uniq (chr1,HS)
mylogit <- glm(tp ~ edit+event+mcs_1+mcs_2+mapq_1+mapq_2, data = data1, family = "binomial")
// Mmaps (chr1,HS)

// Ties (chr1,HS)
mylogit <- glm(tp ~ edit+event+mcs_1+mcs_2+max_region+m2p+mapq_1+mapq_2, data = data1, family = "binomial")




