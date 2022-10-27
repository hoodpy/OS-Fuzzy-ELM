cauchy = function(input, a, alpha){
  output = 1/(1 + alpha * (input - a)^2)
  return(output)
}

cauchy_function = function(instance, a_list, alphas){
  num_features = nrow(alphas)
  num_rules = ncol(alphas)
  cauchy_values = matrix(0, nrow=num_features, ncol=num_rules)
  for(i in 1:num_features){
    for(j in 1:num_rules){
      cauchy_values[i, j] = cauchy(instance[i], a_list[j], alphas[i, j])
    }
  }
  cauchy_values = apply(cauchy_values, 2, min)
  cauchy_values = as.data.frame(matrix(cauchy_values, nrow=1, ncol=num_rules))
  return(cauchy_values/sum(cauchy_values))
}

obtain_H_matrix = function(data_extend, G_data){
  H_matrix = data.frame()
  num_samples = nrow(G_data)
  num_rules = ncol(G_data)
  for(i in 1:num_rules){
    h_matrix = as.data.frame(diag(x=G_data[,i], num_samples, num_samples) %*% as.matrix(data_extend))
    if(i==1){
      H_matrix = h_matrix
    } else {
      H_matrix = cbind(H_matrix, h_matrix)
    }
  }
  return(as.matrix(H_matrix))
}

generate_map_matrix = function(variables, a_list, alphas){
  G_data = data.frame()
  for(i in 1:nrow(variables)){
    G_data = rbind(G_data, cauchy_function(variables[i,], a_list, alphas))
  }
  return(obtain_H_matrix(cbind(1, variables), G_data))
}

obtain_label_matrix = function(labels, num_classes){
  label_matrix = matrix(0, nrow=length(labels), ncol=num_classes)
  for(i in 1:num_classes){
    label_matrix[which(labels==i), i] = 1
  }
  return(label_matrix)
}

normial<-function(x){
  return((2*(x-min(x))/(max(x)-min(x)))-1)
}

obtained_acc_G_mean<-function(x){
  the_sum<-0
  the_G_mean<-1
  for(i in 1:nrow(x)){
    the_sum<-the_sum+x[i,i]
    the_G_mean<-the_G_mean*(x[i,i]/sum(x[i,]))
  }
  the_acc<-the_sum/sum(x)
  the_G_mean<-the_G_mean^(1/nrow(x))
  return(list(the_acc*100,the_G_mean*100))
}

model = function(number_rules, num_classes, train_path, samples_number = 0, test_path = "src", single = FALSE, batch_upper = 100) {
  if (single) {
    total_data = read.table(train_path, header = TRUE, sep = ",", stringsAsFactors = TRUE)
  } else {
    data_train = read.table(train_path, header = TRUE, sep = ",", stringsAsFactors = TRUE)
    data_test = read.table(test_path, header = TRUE, sep = ",", stringsAsFactors = TRUE)
    total_data = rbind(data_train, data_test)
    samples_number = nrow(data_train)
  }
  variables_number = ncol(total_data) - 1
  total_data$label = as.numeric(total_data$label)
  total_data$label = as.factor(total_data$label)
  total_data_normial = as.data.frame(lapply(total_data[, c(1:variables_number)], normial))
  total_data = cbind(total_data_normial, total_data[variables_number + 1])
  data = total_data[c(1:samples_number),]
  testing_data = total_data[-c(1:samples_number),]
  categories = seq(1, by = 1, length = num_classes)
  a_list = runif(number_rules, min=-1, max=1)
  alphas = matrix(runif(number_rules*variables_number, min=0.01, max=9.01), nrow=variables_number, ncol=number_rules)
  start_number = 1
  end_number = number_rules * (variables_number + 1)
  training_data = data[c(start_number:end_number),]
  training_data_variables = as.matrix(training_data[,c(1:variables_number)])
  training_data_labels = obtain_label_matrix(training_data[,variables_number+1], num_classes)
  H = generate_map_matrix(training_data_variables, a_list, alphas)
  K = t(H) %*% H + diag(10^-10, ncol(H), ncol(H))
  Beta = solve(K) %*% t(H)%*%training_data_labels
  while(end_number < samples_number){
    start_number = end_number + 1
    end_number = end_number + sample(1:batch_upper,1,replace=TRUE)
    if(end_number > samples_number){
      end_number = samples_number
    }
    training_data = data[c(start_number:end_number),]
    training_data_variables = as.matrix(training_data[,c(1:variables_number)])
    training_data_labels = obtain_label_matrix(training_data[,variables_number+1], num_classes)
    H = generate_map_matrix(training_data_variables, a_list, alphas)
    K = K + t(H) %*% H
    Beta = Beta + solve(K) %*% t(H) %*% (training_data_labels - H %*% Beta)
  }
  testing_data_variables = as.matrix(testing_data[,c(1:variables_number)])
  H = generate_map_matrix(testing_data_variables, a_list, alphas)
  aim_result = as.data.frame(H %*% Beta)
  aim_result$result = 0
  for(i in 1:nrow(aim_result)){
    aim_result[i,ncol(aim_result)] = which.max(aim_result[i,c(1:(ncol(aim_result)-1))])
  }
  table0<-table(testing_data$label,aim_result$result)
  final_result<-obtained_acc_G_mean(table0)
  Acc<-final_result[[1]]
  Gmean<-final_result[[2]]
  comp<-data.frame(Acc,Gmean)
  names(comp)<-c("Acc","Gmean")
  print(comp)
  saver<-read.table("D:/Documents/program/data.csv",header=TRUE,sep=",")
  saver<-rbind(saver,comp)
  write.csv(saver,"D:/Documents/program/data.csv",row.names=FALSE)
}

for (number in 1:50) {
  model(10, 5, "D:/Documents/program/pageblocks_train.csv", 3000, "D:/Documents/program/pageblocks_test.csv")
}
