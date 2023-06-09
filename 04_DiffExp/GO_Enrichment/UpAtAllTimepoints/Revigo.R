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
revigo.data <- rbind(c("GO:0000184","nuclear-transcribed mRNA catabolic process, nonsense-mediated decay",0.07451642654302854,-5.697886373676468,-5.704766577465544,4.350325563996703,-14,0.7238653099420752,0.33006803),
c("GO:0000966","RNA 5'-end processing",0.10857984323843341,-6.284467120928864,-2.6124074829551325,4.513816672940482,-11.42021640338319,0.6000295050054274,0.34236656),
c("GO:0002181","cytoplasmic translation",0.12046428693071752,-6.840134276508772,0.07560382796438224,4.5589244643947815,-27.2839966563652,0.6393539926750533,0.16724244),
c("GO:0006413","translational initiation",0.5201814019133255,-5.906041879645772,-0.7156339807501588,5.1942117566342985,-16.769551078621728,0.5945069912619841,0.47198297),
c("GO:0006614","SRP-dependent cotranslational protein targeting to membrane",0.13951997624308687,1.1252150494119093,5.611687151403271,4.622700906045888,-15.886056647693163,0.7687496050420651,0.0091555),
c("GO:0006725","cellular aromatic compound metabolic process",26.410310601416047,2.0938915301280407,-6.4515740799345265,6.8998277223761,-13,0.8178724357140909,0.28457639),
c("GO:0010467","gene expression",9.706370753664654,-2.3548677612974416,0.8302859108398974,6.465111183714701,-22.2839966563652,0.810587265951598,0.21063862),
c("GO:0034641","cellular nitrogen compound metabolic process",30.187049102942364,2.9438196003455137,-4.919524053054954,6.957874867913695,-17.468521082957746,0.8020685393996184,0.16383891),
c("GO:0034660","ncRNA metabolic process",3.769763500568378,-3.2148814950790436,-4.056250838656359,6.054368647279385,-30,0.6067115028055005,0.01324497),
c("GO:0042254","ribosome biogenesis",1.6536972228253697,-3.3164362397556078,5.8295231427892675,5.696511029454469,-30,0.5834732740656142,0),
c("GO:0043170","macromolecule metabolic process",38.34861513800723,3.8662554791939376,0.33046154188482696,7.061803881678349,-10.113509274827518,0.8749728242436937,0.19244037),
c("GO:0043603","amide metabolic process",6.518961624921768,5.972783383391009,-3.2119967048693163,6.292232804957447,-9.67778070526608,0.8630939630676597,0.18668265),
c("GO:0046483","heterocycle metabolic process",26.396317359150885,1.2584934592820038,-7.22798277507035,6.899597554781196,-13.10790539730952,0.8178918145463975,0.28451677),
c("GO:0090502","RNA phosphodiester bond hydrolysis, endonucleolytic",0.5997436445999855,-4.73229176891593,-4.4356216851751835,5.2560222219724215,-13.886056647693163,0.6723263700817557,0.41206749),
c("GO:1901360","organic cyclic compound metabolic process",27.27731619259278,4.203186057259491,2.218225479318229,6.9138558498071845,-11.67778070526608,0.8840373190871517,0.09491393));

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
p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(one.data$value), 0) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = log_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) ));
p1 <- p1 + scale_size( range=c(5, 30)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.15, ];
p1 <- p1 + geom_text( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 );
p1 <- p1 + labs (x = "semantic space x", y = "semantic space y");
p1 <- p1 + theme(legend.key = element_blank()) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10)


# --------------------------------------------------------------------------
# Output the plot to screen

#p1;

svg("UpInAll.svg")
print(p1)
dev.off()

# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).

# ggsave("/path_to_your_file/revigo-plot.pdf");
