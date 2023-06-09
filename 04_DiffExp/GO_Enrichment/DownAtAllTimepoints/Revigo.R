# A plotting R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# --------------------------------------------------------------------------
# If you don't have the ggplot2 package installed, uncomment the following line:
# install.packages( "ggplot2" );
library( ggplot2 );

# --------------------------------------------------------------------------
# If you don't have the scales package installed, uncomment the following line:
# install.packages( "scales" );
library( scales );

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","plot_X","plot_Y","log_size","value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0000280","nuclear division",0.16930925072058378,-0.8112083829104746,8.042224831076382,4.706743378506736,-10.537602002101044,0.6312627422559939,0.00920546),
c("GO:0000727","double-strand break repair via break-induced replication",0.012795817208009229,-6.975262017049069,1.8669416984859089,3.5852350633657752,-5.958607314841775,0.7650043399833629,0.40623967),
c("GO:0006099","tricarboxylic acid cycle",0.5408768983178572,4.686590497459715,4.493263296617174,5.21115526206114,-6.4089353929735005,0.944244221853742,0.07610317),
c("GO:0006101","citrate metabolic process",0.011764701186568404,4.702787320600529,0.5997831328307804,3.5487578285737045,-6.4089353929735005,0.9562745912409275,0.05519749),
c("GO:0006259","DNA metabolic process",5.803533415077551,-6.964029065546808,4.29007594216371,6.241746897108298,-7.376750709602099,0.8002436924604996,0.37848828),
c("GO:0006260","DNA replication",1.444034747678592,-6.183227622546943,2.8057102400673952,5.637632802972698,-14.920818753952375,0.6843676090836338,0),
c("GO:0006275","regulation of DNA replication",0.10787136674628214,3.16233922236227,-2.7337930353575928,4.510973731348933,-6.886056647693163,0.9661806239762096,-0),
c("GO:0007049","cell cycle",1.795399173617054,0.13185683190905242,-5.21177012470518,5.732215984404404,-12.008773924307505,0.983296041337164,0.01183773),
c("GO:0022402","cell cycle process",0.9753891897530218,-4.795928681254623,-4.615197933981784,5.4672335778545404,-13.142667503568731,0.38930007312316867,0.01102324),
c("GO:0022616","DNA strand elongation",0.04290773121479571,-5.45377129835678,2.1667507996125663,4.110623375233331,-9.207608310501746,0.7483731322101064,0.44722861),
c("GO:0043570","maintenance of DNA repeat elements",0.04744464170913534,-3.720069103792,4.7425737325039234,4.154271775993095,-7.366531544420414,0.5159623883687267,0.45910622),
c("GO:0070058","tRNA gene clustering",0.0002627682764316946,0.1454562058854279,7.413740938265733,1.9030899869919435,-6.3872161432802645,0.783405080624558,0.33505724),
c("GO:1902981","synthesis of RNA primer involved in mitotic DNA replication",3.6587987857577735E-05,-5.513475063114702,-0.8991914399721637,1.0791812460476249,-7.455931955649724,0.38874152273702056,0.40591833));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$log_size <- as.numeric( as.character(one.data$log_size) );
one.data$value <- as.numeric( as.character(one.data$value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
#head(one.data);


# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = value, size = log_size), alpha = I(0.6) );
#p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(one.data$value), 0) );
# Change limit to -30 so colours sync with upregulated genes' plot
p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c(-30, 0) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = log_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) ));
p1 <- p1 + scale_size( range=c(5, 30)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.15, ];
p1 <- p1 + geom_text( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 );
p1 <- p1 + labs (x = "semantic space x", y = "semantic space y");
p1 <- p1 + theme(legend.key = element_blank()) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);


# --------------------------------------------------------------------------
# Output the plot to screen

#p1;

svg("DownInAll.svg")
print(p1)
dev.off()

# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).

# ggsave("/path_to_your_file/revigo-plot.pdf");
