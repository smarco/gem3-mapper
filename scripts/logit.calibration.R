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
 * SE
 */
// None
mylogit <- glm(tp ~ , data = data1, family = "binomial")

// All
mylogit <- glm(tp ~ edit+sub_edit+event+sub_event+swg+sub_swg+mcs+max_region+fs_matches+sub_matches+sub_cand, data = data1, family = "binomial")

// Uniq (chr1,HS)
mylogit <- glm(tp ~ event+swg+mcs, data = data1, family = "binomial")
mylogit <- glm(tp ~ event+swg+mcs+max_region+sub_cand, data = data1, family = "binomial")

// Mmaps (chr1,HS)
mylogit <- glm(tp ~ event+swg+sub_swg+mcs, data = data1, family = "binomial")
mylogit <- glm(tp ~ event+sub_event+swg+sub_swg+mcs+sub_matches, data = data1, family = "binomial")

// D1-Ties (chr1,HS)
mylogit <- glm(tp ~ edit+sub_edit+event+sub_event+swg+sub_swg+mcs+max_region+fs_matches+sub_matches, data = data1, family = "binomial")
mylogit <- glm(tp ~ edit+sub_edit+event+sub_event+swg+sub_swg+mcs+max_region+fs_matches+sub_matches, data = data1, family = "binomial")

/*
 * PE
 */
mylogit <- glm(tp ~ , data = data1, family = "binomial")

// All
mylogit <- glm(tp ~ edit+sub_edit+event+sub_event+swg+sub_swg+mcs+max_region+fs_matches+sub_matches+sigma+sub_sigma+mapq_e1+mapq_e2, data = data1, family = "binomial")

// Uniq (chr1,HS)
mylogit <- glm(tp ~ edit+event+swg+mcs+sigma+mapq_e1, data = data1, family = "binomial")


// Mmaps (chr1,HS)
mylogit <- glm(tp ~ event+mcs, data = data1, family = "binomial")

// Ties
mylogit <- glm(tp ~ edit+sub_edit+swg+sub_swg+fs_matches, data = data1, family = "binomial")







