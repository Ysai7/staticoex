#' Multivariate Analysis Of Multi-Way Data Matrices
#' @param x set of matrices juxtaposed in vertical order containing the data of the P variables.
#' @param y set of matrices juxtaposed in vertical order containing the data of the Q variables.
#' @param condition A factor object representing the categorical variable that will define the third-path for the analysis.
#' @param biplot Logical value, TRUE by default, that will calculate and graphically represent the compromise and intra-structure considering the HJ-Biplot perspective.
#' @param preprocessing Standardization method to be used. preprocessing="B.Total" performs total standardization, preprocessing="B.Partial" partial standardization, preprocessing="A.Norma" <by default> standardization by standard.
#' @param compx Indicates the component of x to be used for the representation of the compromise and intra-structure plots on the bi-factor plane.
#' @param compy Indicates the component of y to be used for the representation of the compromise and intra-structure plots on the bi-factor plane.
#' @return A series of indices (contributions and representation qualities) and compromise and intra-structure plots
#' @export HJSTATICO
#' @examples
#' result=HJSTATICO(x=env_data,y=spe_data,condition=etiq_data$space,biplot=TRUE,preprocessing = "A.Norma",compx=1,compy=2)
#' result$Compromise
#' result$HJB
#' result$HJBK
#' result$CrossTables;result$X;result$Y
#' @import pracma
#' @import MFAg
#' @import ade4
#' @import methods
#' @import graphics

HJSTATICO=function(x,y,condition,biplot,preprocessing,compx,compy) {
 if(!is.logical(biplot)){
  biplot=TRUE
 }
 if(exists("preprocessing")==FALSE){
  preprocessing="A.Norma"
 }
 if (is(condition,"factor")){
  condition<-as.character(condition)
 }

 tp=preprocessing
 dimx=compx
 dimy=compy
 k=length(unique(condition))##dim(table(condition))
 dimenx=dim(x)
 dimeny=dim(y)
 ##n=dimenx[1]/k #aqui supone que todas las ktablas tienen el mismo n---cambiar
 p=dimenx[2]
 q=dimeny[2]

 prep=staticoex::Prep.Tables(x=x,y=y,condition=condition,preprocessing=tp)
 X=prep$X
 Y=prep$Y
 index=prep$index
 #COINERTIA
 for (i in 1:k) {
  dn=1/index[i,4]
  Dn=diag(dn,nrow=index[i,4], ncol=index[i,4])
  tempx=X[index[i,2]:index[i,3],]
  tempx=as.matrix(tempx)
  tempy=Y[index[i,2]:index[i,3],]
  tempy=as.matrix(tempy)
  tempz=t(tempy)%*%Dn%*%tempx
  if (i==1) {
   Z=tempz
   Z1=tempz}#Z1 esta demas
  else  {
   Z=rbind(Z,tempz)
   Z1=cbind(Z1,tempz) }#Z1 esta demas
 }

 Z #matriz de coinercia que contiene las Zk tablas de coinercia

 #crea matriz index guía para tablas Z - Coinertia
 dimenz=dim(Z)
 indexz=matrix(nrow=k,ncol=5)
 indexz=as.data.frame(indexz)
 colnames(indexz)=c("Tables","Start","Finish","rows","columns")
 temp=dimenz[1]/k
 indexz[1,1]=index[1,1]
 indexz[1,2]=1
 indexz[1,3]=temp
 indexz[1,4]=(indexz[1,3]-indexz[1,2])+1
 indexz[1,5]=dimenz[2]
 for (i in 2:k) {
  indexz[i,1]=index[i,1]
  indexz[i,2]=indexz[i-1,3]+1
  indexz[i,3]=temp*i
  indexz[i,4]=(indexz[i,3]-indexz[i,2])+1
  indexz[i,5]=dimenz[2]
 }
 indexz

 dimenz1=dim(Z1)#parece que esta demas
 indexz1=matrix(nrow=k,ncol=5)
 indexz1=as.data.frame(indexz1)
 colnames(indexz1)=c("Tables","Start","Finish","rows","columns")
 temp=dimenz1[2]/k
 indexz1[1,1]=index[1,1]
 indexz1[1,2]=1
 indexz1[1,3]=temp
 indexz1[1,4]=dimenz1[1]
 indexz1[1,5]=dimenz1[2]/k
 for (i in 2:k) {
  indexz1[i,1]=index[i,1]
  indexz1[i,2]=indexz1[i-1,3]+1
  indexz1[i,3]=temp*i
  indexz1[i,4]=dimenz1[1]
  indexz1[i,5]=dimenz1[2]/k
 }
 indexz1#parece que esta demas

 #INTERSTRUCTURE-PTA
 C=matrix(nrow=k,ncol=k)
 for (i in 1:k) {
  for (j in 1:k) {
   temp1=Z[indexz[i,2]:indexz[i,3],]
   temp2=Z[indexz[j,2]:indexz[j,3],]
   C[i,j]=sum(temp1*temp2)
  }
 }
 C #calcular C como vectorizacion y probar

 #RV calculo
 d=diag(C)
 RV=matrix(nrow=k,ncol=k)
 row.names(RV)=index[,1]
 colnames(RV)=index[,1]
 for (i in 1:k) {
  for (j in 1:k){
   RV[i,j]=C[i,j]/(sqrt(d[i])*sqrt(d[j]))
  }
 }

 myLabels=index[,1]
 print("RV")
 print(round(RV,3))
 l=Rank(RV)
 dvs=svd(RV,nu=min(l,l),nv=min(l,l))
 dvs$d
 dvs$u
 U=dvs$u[,1:2]
 if(min(U[,1])<0){
  U=U*-1}
 GG=U%*%sqrt(diag(dvs$d[1:2]))
 u1=(U[,1])
 one=matrix(1,nrow=length(u1),ncol=1)
 alfa=u1*((sum((u1)*one))^-1) #tables weight
 tb1=matrix(nrow=length(dvs$d),ncol=4) #tabla de eigenvalues, %var y % var acumulada
 colnames(tb1)=c("No.","Eigenvalues","% Inercia", "% Acum Inercia")
 tb1[,2]=round(t(dvs$d),3)
 for (i in 1:length(dvs$d)) {
  tb1[i,1]=i
  tb1[i,3]=round(((dvs$d[i]/sum(dvs$d))*100),3)
  if (i==1){
   tb1[i,4]=round(tb1[i,3],3) }
  else {
   tb1[i,4]=round((tb1[i,3]+tb1[i-1,4]),3) }
 }
 print("Interstructure")
 print(tb1)

 #Obtener el compromiso
 compromise=matrix(nrow=indexz[1,4],ncol=indexz[1,5])
 for(i in 1:k) {
  temp=matrix(nrow=indexz[1,4],ncol=indexz[1,5])
  temp=Z[indexz[i,2]:indexz[i,3],]
  temp=as.matrix(temp)
  temp2=alfa[i]*temp #este si va por Abdi
  if (i==1){
   compromise=temp2 }
  else{
   compromise=compromise+temp2
  }
 }

 #obtener coseno2
 C1=matrix(nrow=k,ncol=1)
 for (i in 1:k) {
  temp1=Z[indexz[i,2]:indexz[i,3],]
  temp2=compromise
  C1[i,1]=sum(temp1*temp2)
 }
 C1
 #coseno2 calculo
 varcomp=sum(compromise*compromise)
 cos2=matrix(nrow=k,ncol=1)
 row.names(cos2)=index[,1]
 for (i in 1:k) {
  cos2[i,1]=C1[i,1]/(sqrt(d[i])*sqrt(varcomp))
 }

 #analisis compromiso con gsvd SIP
 dimenco=dim(compromise)
 csvd=svd(compromise,nu=min(Rank(compromise),Rank(compromise)),nv=min(Rank(compromise),Rank(compromise)))
 delta=diag(csvd$d)
 F=compromise%*%csvd$v #SI coincide (articulo Puri HJ biplot, marcadores filas)
 Q=t(compromise)%*%csvd$u #SI coincide (articulo Puri HJ biplot, marcadores columnas)

 #summary compromise
 delta2=(csvd$d)^2
 tb2=matrix(nrow=length(csvd$d),ncol=4)
 colnames(tb2)=c("No.","Eigenvalues","% Inercia","% Acum Inercia")
 for (i in 1:length(csvd$d)){
  tb2[i,1]=i
  tb2[i,2]=round(delta2[i],3)
  tb2[i,3]=round(((delta2[i]/sum(delta2))*100),3)
  if (i==1){
   tb2[i,4]=round(tb2[i,3],3) }
  else {
   tb2[i,4]=round((tb2[i,3]+tb2[i-1,4]),3) }
 }
 print("Compromise")
 print(tb2)

 tb11=matrix(nrow=k,ncol=3)
 colnames(tb11)=c("Tables","Weights","Cos2")
 for(i in 1:k){
  tb11[i,1]=index[i,1]
  tb11[i,2]=round(alfa[i],3) #lo correcto es alfa[i]
  tb11[i,3]=round(cos2[i,1],3) }
 tb11=as.data.frame(tb11)
 print("K-tables Indexs")
 print(tb11)

 ####Interestructure, Compromise, weigths-cos2 graphics######
 #SALIDAS GRAFICAS RESULTADOS
 plothjstat=function(dimxc=1,dimyc=2){
  mycolorI=as.vector(rep("gray",length(dvs$d)))
  for (i in 1:length(dvs$d)) {
   if (i== 1)   {
    mycolorI[i]="black" }
   if (i == 2){
    mycolorI[i]="black"  }
  }
  mycolorC=as.vector(rep("gray",length(csvd$d)))
  for (i in 1:length(csvd$d)) {
   if (i== dimxc)   {
    mycolorC[i]="black" }
   if (i ==dimyc){
    mycolorC[i]="black"  }
  }
  par(mfrow=c(2,2))
  s.corcircle(cbind(GG[,1],GG[,2]),sub="Interstructure",possub = "topright",csub=1.5,label="")
  text(GG[,1]+0.08, GG[,2], labels=myLabels, cex=0.8,col="black")
  plot(alfa,cos2,xlab="Weights",ylab="Cos2",pch=19,main ="Index",col="black",xlim=c(min(alfa)-0.01,max(alfa)+0.01),ylim=c(min(cos2)-0.01,max(cos2)+0.01)) #no es u1, es alfa
  text(alfa,cos2,labels=index[,1],adj=1)#no es u1, lo correcto es alfa, solo estoy probando
  barplot(tb1[,3],col=mycolorI,main ="Eigenvalues Inter-structure",ylab="%",axes=TRUE,names.arg = tb1[,1])
  barplot(tb2[,3],col=mycolorC,main="Eigenvalues Compromise", ylab="%",axes=TRUE,names.arg = tb2[,1])
 }
 plothjstat(1,2)
 HJB=hjb(biplot=biplot,F=F,Q=Q,dimx=compx,dimy=compy,compromise=compromise)
 graphcomp(dimx=compx,dimy=compy,biplot=biplot,F=F,Q=Q,compromise=compromise,HJB=HJB)


 #codigo interactivo intraestructura - SI
 FK=matrix(nrow=length(Z1[,1]),ncol=k*length(csvd$v[,1]))
 QK=matrix(nrow=k*length(Z1[,1]),ncol=length(csvd$u[,1]))
 for(i in 1:k){
  tempf=matrix(nrow=indexz1[i,4],ncol=length(csvd$v[,1]))
  tempf=(Z1[,indexz1[i,2]:indexz1[i,3]])%*%csvd$v
  tempq=matrix(nrow=indexz1[i,5],ncol=length(csvd$u[,1]))
  tempq=t(Z1[,indexz1[i,2]:indexz1[i,3]])%*%(csvd$u)
  if(i==1) {
   FK=tempf
   QK=tempq }
  else {
   FK=cbind(FK,tempf)
   QK=rbind(QK,tempq) }
 }

 FK
 dim(FK)
 QK
 dim(QK)

 #grafico de las trayectorias-intraestructuras
 intra=graphintra(dimx=compx,dimy=compy,FK=FK,QK=QK,biplot=biplot,indexz1=indexz1,compromise=compromise,k=k)
 if (biplot=="FALSE") {
  Lst=list(X=round(X,3),Y=round(Y,3),CrossTables=round(Z,3),Compromise=round(compromise,3))}
 else {
  Lst=list(X=round(X,3),Y=round(Y,3),CrossTables=round(Z,3),Compromise=round(compromise,3),HJB=HJB,HJBK=intra)
 }

}
hjb=function(biplot,F,Q,dimx,dimy,compromise){
 ii=dimx; jj=dimy
 bp=biplot
 J=F^2
 dimJ=dim(J)
 H=Q^2
 dimH=dim(H)
 rowJ=as.vector(rep(1,dimJ[1]))
 alfaJ=as.vector(rep(1,dimJ[2]))
 rowH=as.vector(rep(1,dimH[1]))
 alfaH=as.vector(rep(1,dimH[2]))
 for (i in 1:dimJ[1]){rowJ[i]=sum(J[i,])}
 for (i in 1:dimJ[2]) {alfaJ[i]=sum(J[,i])}
 for (i in 1:dimH[1]){rowH[i]=sum(H[i,])}
 for (i in 1:dimH[2]) {alfaH[i]=sum(H[,i])}
 CRTi=as.vector(rep(1,dimJ[1]));  #contribucion relativa a la traza del elemento
 CREiFl=matrix(nrow=dimJ[1],ncol=2); row.names(CREiFl)=row.names(compromise); names(CREiFl[1])=ii; names(CREiFl[2])=jj
 CRFlEi=matrix(nrow=dimJ[1],ncol=2); row.names(CRFlEi)=row.names(compromise); names(CRFlEi[1])=ii; names(CRFlEi[2])=jj
 row
 for (i in 1:dimJ[1]) {CRTi[i]=rowJ[i]/sum(alfaJ)}
 for (i in 1:dimJ[1]){
  CREiFl[i,1]=J[i,ii]/alfaJ[ii]; CREiFl[i,2]=J[i,jj]/alfaJ[jj]
  CRFlEi[i,1]=J[i,ii]/rowJ[i]; CRFlEi[i,2]=J[i,jj]/rowJ[i] }
 CRTj=as.vector(rep(1,dimH[1]))#contribucion relativa a la traza de la variable
 CREjFl=matrix(nrow=dimH[1],ncol=2); row.names(CREjFl)=colnames(compromise); names(CREjFl[1])=ii; names(CREjFl[2])=jj
 CRFlEj=matrix(nrow=dimH[1],ncol=2); row.names(CRFlEj)=colnames(compromise); names(CRFlEj[1])=ii; names(CRFlEj[2])=jj
 for (i in 1:dimH[1]) {CRTj[i]=rowH[i]/sum(alfaH)}
 for(i in 1:dimH[1]){
  CREjFl[i,1]=H[i,ii]/alfaH[ii]; CREjFl[i,2]=H[i,jj]/alfaH[jj]
  CRFlEj[i,1]=H[i,ii]/rowH[i];CRFlEj[i,2]=H[i,jj]/rowH[i] }
 QLREi=matrix(nrow=dimJ[1],ncol=1)
 QLREj=matrix(nrow=dimH[1],ncol=1)
 QLREi=CRFlEi[,1]+CRFlEi[,2]; names(QLREi)=row.names(CRFlEi)  # names(QLREi)=c("QLREi")
 QLREj=CRFlEj[,1]+CRFlEj[,2]; names(QLREj)=row.names(CRFlEj) #names(QLREi)=c("QLREj")
 Lst=list(CRTi=round(CRTi,3),CREiFl=round(CREiFl,3),CRFlEi=round(CRFlEi,3),CRTj=round(CRTj,3),CREjFl=round(CREjFl,3),CRFlEj=round(CRFlEj,3),QLREi=round(QLREi,3),QLREj=round(QLREj,3))
}
graphcomp=function (dimx,dimy,F,Q,biplot,compromise,HJB) {
 x=paste("Axis",dimx,sep="")
 y=paste("Axis",dimy,sep="")
 i=dimx; j=dimy
 maxxQ=max(abs(Q[,i])); maxyQ=max(abs(Q[,j]))
 maxxF=max(abs(F[,i])); maxyF=max(abs(F[,j]))
 if (biplot=="FALSE") {
  par(mfrow=c(1,2))
  plot(Q[,i],Q[,j],pch=20,col="white",main="Compromise",xlab=x, ylab=y,xlim=c(-maxxQ-0.05,maxxQ),ylim=c(-maxyQ,maxyQ))
  text(Q[,i],Q[,j],labels=colnames(compromise))
  arrows(x0=0,y0=0,x1=Q[,i],y1=Q[,j],length=0.1)
  abline(h=0,v=0)
  plot(F[,i],F[,j],pch=20,col="white",main="Compromise",xlab=x, ylab=y,xlim=c(-maxxF-0.05,maxxF),ylim=c(-maxyF,maxyF))
  text(F[,i],F[,j],labels=row.names(compromise),adj=1)
  arrows(x0=0,y0=0,x1=F[,i],y1=F[,j],length=0.1)
  abline(h=0,v=0) }
 else {
  par(mfrow=c(1,1))
  qq=cbind(Q[,i],Q[,j])
  ff=cbind(F[,i],F[,j])
  tt=rbind(qq,ff)
  maxb=max(maxxQ,maxxF)
  mayb=max(maxyQ,maxyF)
  plot(qq[,1],qq[,2],pch=20,col="white",main="Compromise:HJBiplot",xlab=x, ylab=y,xlim=c(-maxb-0.05,maxb),ylim=c(-mayb,mayb))
  arrows(x0=0,y0=0,x1=qq[,1],y1=qq[,2],length=0.1,col="blue")
  text(qq[,1],qq[,2],labels=colnames(compromise),col="blue",adj=1)
  points(ff[,1],ff[,2],pch=16,cex=HJB$QLREi,col="black")
  text(ff[,1],ff[,2],labels=row.names(compromise),cex=HJB$QLREi,adj=1)
  abline(h=0,v=0)
 }
}
graphintra=function(dimx,dimy,FK,QK,indexz1,biplot,compromise,k) {
 dimx; dimy
 d1=dimx; d2=dimy
 x=paste("Axis",dimx,sep="")
 y=paste("Axis",dimy,sep="")
 if (biplot=="FALSE") {
  for (i in 1:k) {
   par(mfrow=c(1,2))
   tempFK=FK[,indexz1[i,2]:indexz1[i,3]]
   tempQK=QK[indexz1[i,2]:indexz1[i,3],]
   maxxQK=max(abs(tempQK[,d1])); maxyQK=max(abs(tempQK[,d2]))
   maxxFK=max(abs(tempFK[,d1])); maxyFK=max(abs(tempFK[,d2]))
   maxx=max(maxxQK,maxxFK); mayy=max(maxyQK,maxyFK)
   plot(tempQK[,d1],tempQK[,d2],pch=20,col="white",xlab=x,ylab=y,xlim = c(-maxxQK-0.07,maxxQK),ylim=c(-maxyQK,maxyQK),main=indexz1[i,1])
   arrows(x0=0,y0=0,x1=tempQK[,d1],y1=tempQK[,d2],length=0.1)
   text(tempQK[,d1],tempQK[,d2],labels=colnames(compromise),adj=1)
   abline(h=0,v=0)
   plot(tempFK[,d1],tempFK[,d2],pch=20,col="white",xlab=x,ylab=y,xlim = c(-maxxFK-0.07,maxxFK),ylim=c(-maxyFK,maxyFK),main=indexz1[i,1])
   arrows(x0=0,y0=0,x1=tempFK[,d1],y1=tempFK[,d2],length=0.1)
   text(tempFK[,d1],tempFK[,d2],labels=row.names(compromise),adj=1)
   abline(h=0,v=0) }
 }
 else {
  for(i in 1:k) {
   par(mfrow=c(1,1))
   tempFK=FK[,indexz1[i,2]:indexz1[i,3]]; F=tempFK
   tempQK=QK[indexz1[i,2]:indexz1[i,3],]; Q=tempQK
   HJBK=hjb(biplot=TRUE, F=F, Q=Q,dimx=d1,dimy=d2,compromise=compromise)
   qq=cbind(tempQK[,d1],tempQK[,d2])
   ff=cbind(tempFK[,d1],tempFK[,d2])
   tt=rbind(qq,ff)
   maxxQK=max(abs(tempQK[,d1])); maxyQK=max(abs(tempQK[,d2]))
   maxxFK=max(abs(tempFK[,d1])); maxyFK=max(abs(tempFK[,d2]))
   maxx=max(maxxQK,maxxFK); mayy=max(maxyQK,maxyFK)
   plot(qq[,1],qq[,2],pch=20,col="white",main=paste(indexz1[i,1],"HJBiplot",sep=":"),xlab=x,ylab=y,xlim=c(-maxx-0.05,maxx),ylim=c(-mayy,mayy))
   arrows(x0=0,y0=0,x1=qq[,1],y1=qq[,2],length=0.1,col="blue")
   text(qq[,1],qq[,2],labels=colnames(compromise),col="blue",adj=1)
   points(ff[,1],ff[,2],pch=16,cex=HJBK$QLREi,col="black")
   text(ff[,1],ff[,2],labels=row.names(compromise),cex=HJBK$QLREi,adj=1)
   abline(h=0,v=0)
   if(i==1){
    t1=0;t2=0;t3=0;t4=0;t5=0;t6=0;t7=0;t8=0
    t1=HJBK$CRTi; t2=HJBK$CREiFl; t3=HJBK$CRFlEi;
    t4=HJBK$CRTj; t5=HJBK$CREjFl; t6=HJBK$CRFlEj;
    t7=HJBK$QLREi; t8=HJBK$QLREj }
   else{
    t1=c(t1,HJBK$CRTi); t2=rbind(t2,HJBK$CREiFl); t3=rbind(t3,HJBK$CRFlEi);
    t4=c(t4,HJBK$CRTj); t5=rbind(t5,HJBK$CREjFl); t6=rbind(t6,HJBK$CRFlEj);
    t7=c(t7,HJBK$QLREi); t8=c(t8,HJBK$QLREj) }}
  Lst=list(CRTi=round(t1,3),CREiFl=round(t2,3),CRFlEi=round(t3,3),CRTj=round(t4,3),CREjFl=round(t5,3),CRFlEj=round(t6,3),QLREi=round(t7,3),QLREj=round(t8,3))
 }
}
