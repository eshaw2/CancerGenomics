ref_overlap_ct = function(data, all_TPs, all_FPs, our_full_cosmic, label=F){
  our_TPs = all_TPs %>% mutate(ref_grp = 'TP')
  our_FPs = all_FPs %>% mutate(ref_grp = 'FP')
  all_drivers = bind_rows(our_TPs,our_FPs) %>% select(Codon_Id,comparison) %>% unique()
  our_TNs = anti_join(data,all_drivers, by='Codon_Id') %>% anti_join(our_full_cosmic, by='Codon_Id')  %>% mutate(ref_grp = 'TN',comparison='UNK')
  
  if (label) {
    return(
      bind_rows(our_TPs%>% mutate(ref_label = "positives"),
                bind_rows(our_FPs,our_TNs) %>% mutate(ref_label = "negatives")
                ) %>%
        group_by(ref_label) %>%
        tally(name="ref_records")
    )
  }
  return( bind_rows(our_TPs,our_FPs,our_TNs) %>% 
            group_by(ref_grp,comparison) %>% 
            tally(name="ref_records"))
}

compare_refs = function(df,all_TPs = read.csv('data/all_TPs.csv'),
                        all_FPs = read.csv('data/potential_FP.csv'), 
                        our_full_cosmic = read.csv('data/full_cosmic.csv')) {
    our_cBios = all_TPs %>% filter(comparison=='chang_TP')
    our_changs = all_TPs %>% filter(comparison=='chang_pred')
    our_browns = all_TPs %>% filter(comparison=='brown_TP')
    our_mutagenes_TP = all_TPs %>% filter(comparison=='li_TP')
    consensus_TPs = inner_join(
                      inner_join(
                        inner_join(our_cBios,our_changs, by='Codon_Id'),
                        our_browns, by='Codon_Id'),
                      our_mutagenes_TP, by='Codon_Id') %>% mutate(comparison='consensus')
    our_mutagenes_FP = all_FPs %>% filter(comparison=='li_FP')
    our_brown_FP = all_FPs %>% filter(comparison=='brown_FP')

    cBio_results = df %>% semi_join(our_cBios, by="Codon_Id") %>% mutate(comparison='chang_TP')
    chang_results = df %>% semi_join(our_changs, by="Codon_Id") %>% mutate(comparison='chang_pred')
    brown_TP_results = df %>% semi_join(our_browns, by="Codon_Id") %>% mutate(comparison='brown_TP')
    brown_FP_results = df %>% semi_join(our_brown_FP, by="Codon_Id") %>% mutate(comparison='brown_FP')
    mutagene_TP_results = df %>% semi_join(our_mutagenes_TP, by="Codon_Id") %>% mutate(comparison='li_TP')
    mutagene_FP_results = df %>% semi_join(our_mutagenes_FP, by="Codon_Id") %>% mutate(comparison='li_FP')
    cosmic_results = df %>% anti_join(our_full_cosmic, by="Codon_Id") %>% mutate(comparison='UNK')
    consensus_results = df %>% semi_join(consensus_TPs, by='Codon_Id')  %>% mutate(comparison='consensus')
    
    comb_results = bind_rows(cBio_results,chang_results,
                             brown_TP_results,brown_FP_results,
                             mutagene_TP_results,mutagene_FP_results) 
    
    TP_results = bind_rows(cBio_results, mutagene_TP_results, consensus_results,
                           chang_results, brown_TP_results) %>% 
                  mutate(ref_grp = 'TP')
    
    FP_results = bind_rows(mutagene_FP_results,brown_FP_results) %>% 
                    mutate(ref_grp = 'FP')
    
    TN_results = anti_join(df,comb_results, by='Codon_Id') %>%
                    semi_join(cosmic_results, by='Codon_Id') %>%
                    mutate(comparison='UNK',
                           ref_grp = 'TN')
    all_results = bind_rows(TP_results, FP_results, TN_results) %>%
                    mutate(comparison = factor(comparison),
                           ref_grp = factor(ref_grp,
                                            levels=c("FP","TP","TN"))
                    )
    div_results = bind_rows( 
                        TP_results%>% mutate(ref_label = "positives"),
                        bind_rows(TN_results,FP_results) %>% mutate(ref_label = "negatives")
                      )
    
    return(list(TP = TP_results, FP = FP_results, 
                TN = TN_results, 
                all = all_results, div = div_results))
}