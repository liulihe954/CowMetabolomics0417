
str(step5_2)

step5_2_test = step5_2
step5_2_test$plot_order_x = 5

Final_edgelist_low_test = step5_2_test


#Final_edgelist_high = step6_2


write.csv(Final_edgelist_low_test,"Final_edgelist_cyto_low_test_attri.txt",quote = F)
# write.csv(Final_edgelist_high,"Final_edgelist_cyto_high.txt",quote = F)
