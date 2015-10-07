package sa;

/**
 *  Pairwise Optimal Sequence Alignment
 *
 *  This program uses the dynamic programming based sequence alignment algorithm, 
 *   with the following three variations: global alignment, local alignment, and global alignment 
 *   with affine gap penalties.
 *
 *  This takes O(n^2) time to construct the matrix and 
 *   linear time to extract the optimal alignment from the matrix constructed. 
 *   The alignment problems are only for DNA sequences.
 *  
 *  @version 1.0 03 March 2014
 *  @author Cedric Jo
 */
public class SeqAlignMethod {
    
    /** String for global alignment output */
    public StringBuilder g;
    
    /** String for local alignment output */
    public StringBuilder l;
    
    /** String for global alignment with affine penalty output */
    public StringBuilder af;
    
    
    /** 
     *  Global sequence alignment function
     *  
     *    Returns global aligned DNA sequences (string)
     *    Parameters:
     *      s1 = score for perfect match
     *      s2 = score for same base substitution
     *      s3 = score for other substitution
     *      p1 = gap penalty
     *      sq1, sq2 = DNA sequences
     */
    public String global(int s1, int s2, int s3, int p1, String sq1, String sq2) {
        
        /** matrix for dynamic programming */
        int[][] gOpt = new int[sq2.length()+1][sq1.length()+1];
        
        /** gap penalty */
        int gap = -1*Math.abs(p1);
       
        Integer a = null;
        
        
        /** Filling first row & column with gap penalty */
        for (int i=1; i <= sq2.length(); i++) {
            gOpt[i][0] = gOpt[i-1][0] + gap;
        }
        for (int j=1; j <= sq1.length(); j++) {
            gOpt[0][j] = gOpt[0][j-1] + gap;
        }
        
        /** Filling the rest of the matrix */
        for (int i=1; i <= sq2.length(); i++) {
            for (int j=1; j <= sq1.length(); j++) {
                
                /** diagonal score */
                int scoreDiag;
                
                if (sq2.charAt(i-1) == sq1.charAt(j-1)) {
                    scoreDiag = gOpt[i-1][j-1] + s1;
                }
                else if (sq2.charAt(i-1) == 'A' && sq1.charAt(j-1) == 'G') {
                    scoreDiag = gOpt[i-1][j-1] + s2;
                }
                else if (sq2.charAt(i-1) == 'G' && sq1.charAt(j-1) == 'A') {
                    scoreDiag = gOpt[i-1][j-1] + s2;
                }
                else if (sq2.charAt(i-1) == 'C' && sq1.charAt(j-1) == 'T') {
                    scoreDiag = gOpt[i-1][j-1] + s2;
                }
                else if (sq2.charAt(i-1) == 'T' && sq1.charAt(j-1) == 'C') {
                    scoreDiag = gOpt[i-1][j-1] + s2;
                }
                else {
                    scoreDiag = gOpt[i-1][j-1] + s3;
                }
                
                int scoreLeft = gOpt[i][j-1] + gap;
                int scoreUp = gOpt[i-1][j] + gap;
                
                gOpt[i][j] = Math.max(Math.max(scoreDiag, scoreLeft), scoreUp);
            }
        }
        
        /*-----------------------------------------------------------------*/
        /** Print Global alignment table */  
          g = new StringBuilder("");
          g.append("\t" + "  -" + "\t");
          for (int i=0; i <= sq2.length(); i++) 
          {
              if (i==0) {
                for (int j=0; j < sq1.length(); j++)
                {
                    g.append("  " + sq1.charAt(j) + "\t");
                }
                g.append("\n");}
              
              
              if (i==0) {
                      g.append("-" + "\t");
                  }
              for (int j=0; j <= sq1.length(); j++)
              {
                  if(i >= 1 && j==0)
                  {
                      g.append(sq2.charAt(i-1) + "\t");
                  }
                  
                  
                  if (gOpt[i][j] >= 0)
                  {
                      g.append("  " + gOpt[i][j] + "\t");
                  }
                  else {
                  g.append(gOpt[i][j] + "\t");
                  }         
                  //System.out.print(gOpt[i][j]+ "\t");
              }
              g.append("\n");
              //System.out.println();
          }  
        return Traceback(gOpt, sq1, sq2, gap, s1, s2, s3);
    }
    
    
    
    
    /**  
     *  Local sequence alignment
     *    
     *    Returns aligned DNA substrings
     *    parameters:
     *      s1 = score for perfect match
     *      s2 = score for same base substitution
     *      s3 = score for other substitution
     *      p1 = gap penalty
     *      sq1, sq2 = DNA sequences
     */
    public String local(int s1, int s2, int s3, int p1, String sq1, String sq2) {
        int[][] lOpt = new int[sq2.length()+1][sq1.length()+1];
        int gap = -1*Math.abs(p1);
        
       
        /** Filling first row & column with 0 */
        for (int i=1; i <= sq2.length(); i++) {
            lOpt[i][0] = 0;
        }
        for (int j=1; j <= sq1.length(); j++) {
            lOpt[0][j] = 0;
        }
        
        /** Filling rest of the matrix */
        for (int i=1; i <= sq2.length(); i++) {
            for (int j=1; j <= sq1.length(); j++) {
                int scoreDiag;
                
                if (sq2.charAt(i-1) == sq1.charAt(j-1)) {
                    scoreDiag = lOpt[i-1][j-1] + s1;
                }
                else if (sq2.charAt(i-1) == 'A' && sq1.charAt(j-1) == 'G') {
                    scoreDiag = lOpt[i-1][j-1] + s2;
                }
                else if (sq2.charAt(i-1) == 'G' && sq1.charAt(j-1) == 'A') {
                    scoreDiag = lOpt[i-1][j-1] + s2;
                }
                else if (sq2.charAt(i-1) == 'C' && sq1.charAt(j-1) == 'T') {
                    scoreDiag = lOpt[i-1][j-1] + s2;
                }
                else if (sq2.charAt(i-1) == 'T' && sq1.charAt(j-1) == 'C') {
                    scoreDiag = lOpt[i-1][j-1] + s2;
                }
                else {
                    scoreDiag = lOpt[i-1][j-1] + s3;
                }
                
                int scoreLeft = lOpt[i][j-1] + gap;
                int scoreUp = lOpt[i-1][j] + gap;
                
                int max = Math.max(Math.max(scoreDiag, scoreLeft), scoreUp);
                    if (max < 0) {
                       lOpt[i][j] = 0;
                    }
                    else {
                       lOpt[i][j] = max;
                    }
            }
        }
        
          // Print local alignment score table
          l = new StringBuilder("");
          l.append("\t" + "  -" + "\t");
          for (int i=0; i <= sq2.length(); i++) {
              if (i==0) {
                for (int j=0; j < sq1.length(); j++) {
                    l.append("  " + sq1.charAt(j) + "\t");
                }
                l.append("\n");}
              
              
              if (i==0) {
                      l.append("-" + "\t");
                  }
              for (int j=0; j <= sq1.length(); j++) {
                  if(i >= 1 && j==0) {
                      l.append(sq2.charAt(i-1) + "\t");
                  }
                  
                  
                  if (lOpt[i][j] >= 0) {
                      l.append("  " + lOpt[i][j] + "\t");
                  }
                  else {
                  l.append(lOpt[i][j] + "\t");
                  }         
              }
              l.append("\n");
          } 
          
        /*-----------------------------------------------------*/
        /** Finding optimal score for local alignment */
        
        /** Find location of max score */
        int tmpI = 0;
        int tmpJ = 0;
        int max = lOpt[0][0];
        for (int i=0; i <= sq2.length(); i++) {
            int[] inner = lOpt[i];
            for (int j=0; j <= sq1.length(); j++) {
                if(inner[j] > max) {
                    max = inner[j];
                    tmpI = i; tmpJ = j;
                }
            }
        }

        /* --------------------------------------------------------------------*/
        // TraceBack from the max score
        String finals1 = "";
        String finals2 = "";
        String s1a = "";
        String s2a = "";
        char ch1 = ' ';
        char ch2 = ' ';
        
        for (int i=tmpI, j=tmpJ; i >= 0 || j >= 0;) {
            int cur = lOpt[i][j];
   
            if (cur==0) {
                break;
            }
            
            else {
            int left = lOpt[i][j-1];
            int dia = lOpt[i-1][j-1];
            int up = lOpt[i-1][j];

            /** Diagonal */
            if (sq2.charAt(i-1) == sq1.charAt(j-1)) {
                ch1 = sq1.charAt(j-1);
                ch2 = sq2.charAt(i-1);

                s1a += ch1;
                s2a += ch2;
                i--;
                j--;
            }
            else if ((sq2.charAt(i-1) == 'A' && sq1.charAt(j-1) == 'G') && dia==cur-s2) {
                ch1 = sq1.charAt(j-1);
                ch2 = sq2.charAt(i-1);

                s1a += ch1;
                s2a += ch2;
                i--;
                j--;
            }
            else if ((sq2.charAt(i-1) == 'G' && sq1.charAt(j-1) == 'A') && dia==cur-s2) {
                ch1 = sq1.charAt(j-1);
                ch2 = sq2.charAt(i-1);

                s1a += ch1;
                s2a += ch2;
                i--;
                j--;
            }
            else if ((sq2.charAt(i-1) == 'C' && sq1.charAt(j-1) == 'T') && dia==cur-s2) {
                ch1 = sq1.charAt(j-1);
                ch2 = sq2.charAt(i-1);

                s1a += ch1;
                s2a += ch2;
                i--;
                j--;
            }
            else if ((sq2.charAt(i-1) == 'T' && sq1.charAt(j-1) == 'C') && dia==cur-s2) {
                ch1 = sq1.charAt(j-1);
                ch2 = sq2.charAt(i-1);

                s1a += ch1;
                s2a += ch2;
                i--;
                j--;
            }
            else if (cur-s3 == dia) {
                ch1 = sq1.charAt(j-1);
                ch2 = sq2.charAt(i-1);

                s1a += ch1;
                s2a += ch2;
                i--;
                j--;
            }
            
            /** Up */
            else if (cur-gap == up) {
                ch2 = sq2.charAt(i-1);
                ch1 = '-';

                s1a += ch1;
                s2a += ch2;
                i--;
            }
                       
            /** Left */
            else // (cur-p1 == left)
            {
                ch1 = sq1.charAt(j-1); 
                ch2 = '-';

                s1a += ch1;
                s2a += ch2;
                j--;
                } 
            }}
        
        finals1 = new StringBuffer(s1a).reverse().toString();
        finals2 = new StringBuffer(s2a).reverse().toString();
        
        String result = finals1 + "\n" + finals2;
        return result;
    } 

    
    
    
    /**  
     *  Affine global sequence alignment
     *    Returns aligned DNA sequences
     *  
     *    Parameters:
     *      s1 = score for perfect match
     *      s2 = score for same base substitution
     *      s3 = score for other substitution
     *      p1 = gap penalty
     *      p2 = affine gap penalty (introduction)
     *      p3 = affine gap penalty (additional)
     *      sq1, sq2 = DNA sequences
     */
    public String affine(int s1, int s2, int s3, int p2, int p3, String sq1, String sq2) {
        int[][] afOpt = new int[sq2.length()+1][sq1.length()+1];    // max score table
        int[][] diaOpt = new int[sq2.length()+1][sq1.length()+1];   // diagonal score table
        int[][] leftOpt = new int[sq2.length()+1][sq1.length()+1];  // (from)left score table
        int[][] upOpt = new int[sq2.length()+1][sq1.length()+1];    // (from)up score table
        
        int a = -1*Math.abs(p2);
        int b = -1*Math.abs(p3);
        
        /*----- Filling first row & column of the tables -----*/
        afOpt[1][0] = a;
        afOpt[0][1] = a;
        leftOpt[1][0] = a;
        upOpt[0][1] = a;
        
        for (int i=2; i <= sq2.length(); i++) 
        {
            afOpt[i][0] = afOpt[i-1][0] + b;
            leftOpt[i][0] = afOpt[i-1][0] + b;
        }
        for (int j=2; j <= sq1.length(); j++)
        {
            afOpt[0][j] = afOpt[0][j-1] + b;
            upOpt[0][j] = afOpt[0][j-1] + b;
        }
        
        /*----- Filling rest of the tables -----*/
        for (int i=1; i <= sq2.length(); i++)
        {
            for (int j=1; j <= sq1.length(); j++)
            {
                int scoreDiag;
                int scoreLeft;
                int scoreUp;
                
                if (sq2.charAt(i-1) == sq1.charAt(j-1))
                {
                    scoreDiag = afOpt[i-1][j-1] + s1;
                }
                else if (sq2.charAt(i-1) == 'A' && sq1.charAt(j-1) == 'G')
                {
                    scoreDiag = afOpt[i-1][j-1] + s2;
                }
                else if (sq2.charAt(i-1) == 'G' && sq1.charAt(j-1) == 'A')
                {
                    scoreDiag = afOpt[i-1][j-1] + s2;
                }
                else if (sq2.charAt(i-1) == 'C' && sq1.charAt(j-1) == 'T')
                {
                    scoreDiag = afOpt[i-1][j-1] + s2;
                }
                else if (sq2.charAt(i-1) == 'T' && sq1.charAt(j-1) == 'C')
                {
                    scoreDiag = afOpt[i-1][j-1] + s2;
                }
                else
                {
                    scoreDiag = afOpt[i-1][j-1] + s3;
                }
                diaOpt[i][j] = scoreDiag;
                
                // Score from left
                if (afOpt[i][j-1] == leftOpt[i][j-1])
                {
                    scoreLeft = afOpt[i][j-1] + b;
                }
                else 
                {
                    scoreLeft = afOpt[i][j-1] + a;
                }
                leftOpt[i][j] = scoreLeft;
                
                // Score from Up
                if(afOpt[i-1][j] == upOpt[i-1][j])
                {
                    scoreUp = afOpt[i-1][j] + b;
                }
                else
                {
                    scoreUp = afOpt[i-1][j] + a;
                }
                upOpt[i][j] = scoreUp; 
                
                afOpt[i][j] = Math.max(Math.max(scoreDiag, scoreLeft), scoreUp);
            }
        }
        
      
         // Print affine score table;
          af = new StringBuilder("");
          af.append("\t" + "  -" + "\t");
          for (int i=0; i <= sq2.length(); i++) 
          {
              if (i==0) {
                for (int j=0; j < sq1.length(); j++)
                {
                    af.append("  " + sq1.charAt(j) + "\t");
                }
                af.append("\n");}
              
              
              if (i==0) {
                      af.append("-" + "\t");
                  }
              for (int j=0; j <= sq1.length(); j++)
              {
                  if(i >= 1 && j==0)
                  {
                      af.append(sq2.charAt(i-1) + "\t");
                  }
                  
                  
                  if (afOpt[i][j] >= 0)
                  {
                      af.append("  " + afOpt[i][j] + "\t");
                  }
                  else {
                  af.append(afOpt[i][j] + "\t");
                  }         
              }
              af.append("\n");
          }
        
        return TracebackA(afOpt, diaOpt, leftOpt, upOpt, sq1, sq2, a, b, s1, s2, s3);
    }    
        
    
    
    
   
    
    
    /**  
     *  Traceback global alignment
     *    Returns aligned DNA sequences
     *  
     *    parameters:
     *      a = score table matrix
     *      sq1, sq2 = DNA sequences 
     */
    public String Traceback(int[][] a, String sq1, String sq2, int p1, int s1, int s2, int s3)
    {
        String finals1 = "";
        String finals2 = "";
        String s1a = ""; 
        String s2a = ""; 
        char ch1 = ' ';
        char ch2 = ' ';
        
        
        for (int i=sq2.length(), j=sq1.length(); i > 0 || j > 0;)
        {
            int cur = a[i][j];
            int left;
            int dia;
            int up;
    
            if (i==0)
            {
                ch1 = sq1.charAt(j-1);           
                ch2 = '-';
                
                s1a += ch1;
                s2a += ch2;
                j--;
            }
            else if(j==0)
            {
               ch1 = '-';     
               ch2 = sq1.charAt(i-1);
                
               s1a += ch1;
               s2a += ch2;
               i--;       
            }
            
            else 
            {
                left = a[i][j-1];
                dia = a[i-1][j-1];
                up = a[i-1][j]; 
          
                // Up
                if (cur-p1 == up)
                {
                    ch2 = sq2.charAt(i-1);
                    ch1 = '-';
                    
                    s1a += ch1;
                    s2a += ch2;
                    i--;
                }

                // Diagonal
                else if (sq2.charAt(i-1) == sq1.charAt(j-1))
                {
                    ch1 = sq1.charAt(j-1);
                    ch2 = sq2.charAt(i-1);
                    
                    s1a += ch1;
                    s2a += ch2;
                    i--;
                    j--;
                }
                else if ((sq2.charAt(i-1) == 'A' && sq1.charAt(j-1) == 'G') && dia==cur-s2)
                {
                    ch1 = sq1.charAt(j-1);
                    ch2 = sq2.charAt(i-1);
                    
                    s1a += ch1;
                    s2a += ch2;
                    i--;
                    j--;
                }
                else if ((sq2.charAt(i-1) == 'G' && sq1.charAt(j-1) == 'A') && dia==cur-s2)
                {
                    ch1 = sq1.charAt(j-1);
                    ch2 = sq2.charAt(i-1);
                    
                    s1a += ch1;
                    s2a += ch2;
                    i--;
                    j--;
                }
                else if ((sq2.charAt(i-1) == 'C' && sq1.charAt(j-1) == 'T') && dia==cur-s2)
                {
                    ch1 = sq1.charAt(j-1);
                    ch2 = sq2.charAt(i-1);
                    
                    s1a += ch1;
                    s2a += ch2;
                    i--;
                    j--;
                }
                else if ((sq2.charAt(i-1) == 'T' && sq1.charAt(j-1) == 'C') && dia==cur-s2)
                {
                    ch1 = sq1.charAt(j-1);
                    ch2 = sq2.charAt(i-1);
                    
                    s1a += ch1;
                    s2a += ch2;
                    i--;
                    j--;
                }
                else if (cur-s3 == dia)
                {
                    ch1 = sq1.charAt(j-1);
                    ch2 = sq2.charAt(i-1);
                    
                    s1a += ch1;
                    s2a += ch2;
                    i--;
                    j--;
                }

                // Left
                else // (cur-p1 == left)
                {
                    ch1 = sq1.charAt(j-1); 
                    ch2 = '-';
                    
                    s1a += ch1;
                    s2a += ch2;
                    j--;
                } 
            }
        }
        finals1 = new StringBuffer(s1a).reverse().toString();
        finals2 = new StringBuffer(s2a).reverse().toString();
        
        String result = finals1 + "\n" + finals2;
        return result; 
    }
    
    
    
    
    /**  
     *  Traceback affine global alignment
     *    
     *    Returns aligned DNA sequences
     *    Parameters:
     *      a = score table maxtrix
     *      sq1, sq2 = DNA sequences 
     */
    public String TracebackA(int[][] affT, int[][] diaT, int[][] leftT, int[][] upT, String sq1, String sq2, int a, int b, int s1, int s2, int s3)
    {
        String finals1 = "";
        String finals2 = "";
        String s1a = ""; 
        String s2a = ""; 
        char ch1 = ' ';
        char ch2 = ' ';
        
        
        for (int i=sq2.length(), j=sq1.length(); i > 0 || j > 0;)
        {
            int cur = affT[i][j];
            int left;
            int dia;
            int up;
    
            if (i==0)
            {
                ch1 = sq1.charAt(j-1);           
                ch2 = '-';
                
                s1a += ch1;
                s2a += ch2;
                j--;
            }
            else if(j==0)
            {
               ch1 = '-';     
               ch2 = sq1.charAt(i-1);
                
               s1a += ch1;
               s2a += ch2;
               i--;       
            }
            
            else 
            {
                left = affT[i][j-1];
                dia = affT[i-1][j-1];
                up = affT[i-1][j]; 
          
                
                // Diagonal
                if (sq2.charAt(i-1) == sq1.charAt(j-1))
                {
                    ch1 = sq1.charAt(j-1);
                    ch2 = sq2.charAt(i-1);
                    
                    s1a += ch1;
                    s2a += ch2;
                    i--;
                    j--;
                }
                else if ((sq2.charAt(i-1) == 'A' && sq1.charAt(j-1) == 'G') && dia==cur-s2)
                {
                    ch1 = sq1.charAt(j-1);
                    ch2 = sq2.charAt(i-1);
                    
                    s1a += ch1;
                    s2a += ch2;
                    i--;
                    j--;
                }
                else if ((sq2.charAt(i-1) == 'G' && sq1.charAt(j-1) == 'A') && dia==cur-s2)
                {
                    ch1 = sq1.charAt(j-1);
                    ch2 = sq2.charAt(i-1);
                    
                    s1a += ch1;
                    s2a += ch2;
                    i--;
                    j--;
                }
                else if ((sq2.charAt(i-1) == 'C' && sq1.charAt(j-1) == 'T') && dia==cur-s2)
                {
                    ch1 = sq1.charAt(j-1);
                    ch2 = sq2.charAt(i-1);
                    
                    s1a += ch1;
                    s2a += ch2;
                    i--;
                    j--;
                }
                else if ((sq2.charAt(i-1) == 'T' && sq1.charAt(j-1) == 'C') && dia==cur-s2)
                {
                    ch1 = sq1.charAt(j-1);
                    ch2 = sq2.charAt(i-1);
                    
                    s1a += ch1;
                    s2a += ch2;
                    i--;
                    j--;
                }
                else if (cur-s3 == dia)
                {
                    ch1 = sq1.charAt(j-1);
                    ch2 = sq2.charAt(i-1);
                    
                    s1a += ch1;
                    s2a += ch2;
                    i--;
                    j--;
                }
                
    
                // Up
                else if (cur-a == up)
                {
                    ch2 = sq2.charAt(i-1);
                    ch1 = '-';
                    
                    s1a += ch1;
                    s2a += ch2;
                    i--;
                }

                // Left
                else if (cur-a == left)
                {
                    ch1 = sq1.charAt(j-1); 
                    ch2 = '-';
                    
                    s1a += ch1;
                    s2a += ch2;
                    j--;
                }
                
                else if (cur-b == up)
                {
                    ch2 = sq2.charAt(i-1);
                    ch1 = '-';
                    
                    s1a += ch1;
                    s2a += ch2;
                    i--;
                }
                else
                {
                    ch1 = sq1.charAt(j-1); 
                    ch2 = '-';
                    
                    s1a += ch1;
                    s2a += ch2;
                    j--;
                }
            }
        }
        finals1 = new StringBuffer(s1a).reverse().toString();
        finals2 = new StringBuffer(s2a).reverse().toString();
        
        String result = finals1 + "\n" + finals2;
        return result; 
    }
}
