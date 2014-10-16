/*  File: SAUI.java
 *  Name: Cedric Jo
 *  Date: 03/03/2014
 * 
 */

package sa;

import javax.swing.*;


public class SAUI extends javax.swing.JFrame {
    SAmethod method = new SAmethod();
    private boolean global;
    private boolean local;
    private boolean affine;
    
    
    public SAUI() {
        initComponents();
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        alignjPanel = new javax.swing.JPanel();
        globaljRadioButton = new javax.swing.JRadioButton();
        localjRadioButton = new javax.swing.JRadioButton();
        affinejRadioButton = new javax.swing.JRadioButton();
        dnasqjPanel = new javax.swing.JPanel();
        dna1jLabel = new javax.swing.JLabel();
        dna2jLabel = new javax.swing.JLabel();
        dna1jTextField = new javax.swing.JTextField();
        dna2jTextField = new javax.swing.JTextField();
        resultjPanel = new javax.swing.JPanel();
        jScrollPane2 = new javax.swing.JScrollPane();
        resultjTextArea = new javax.swing.JTextArea();
        penaltyjPanel = new javax.swing.JPanel();
        p1jLabel = new javax.swing.JLabel();
        p2jLabel = new javax.swing.JLabel();
        p3jLabel = new javax.swing.JLabel();
        p1jTextField = new javax.swing.JTextField();
        p2jTextField = new javax.swing.JTextField();
        p3jTextField = new javax.swing.JTextField();
        scorejPanel = new javax.swing.JPanel();
        s1jTextField = new javax.swing.JTextField();
        s2jTextField = new javax.swing.JTextField();
        s3jTextField = new javax.swing.JTextField();
        s1jLabel = new javax.swing.JLabel();
        s2jLabel = new javax.swing.JLabel();
        s3jLabel = new javax.swing.JLabel();
        excutejPanel = new javax.swing.JPanel();
        runjButton = new javax.swing.JButton();
        clearjButton = new javax.swing.JButton();
        exitjButton = new javax.swing.JButton();
        tablejPanel = new javax.swing.JPanel();
        jScrollPane3 = new javax.swing.JScrollPane();
        tablejTextArea = new javax.swing.JTextArea();

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);

        alignjPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Alignment Type"));

        globaljRadioButton.setText("Global Alignment");
        globaljRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                globaljRadioButtonActionPerformed(evt);
            }
        });

        localjRadioButton.setText("Local Alignment");
        localjRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                localjRadioButtonActionPerformed(evt);
            }
        });

        affinejRadioButton.setText("Affine Global Alignment");
        affinejRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                affinejRadioButtonActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout alignjPanelLayout = new javax.swing.GroupLayout(alignjPanel);
        alignjPanel.setLayout(alignjPanelLayout);
        alignjPanelLayout.setHorizontalGroup(
            alignjPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(alignjPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(alignjPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(globaljRadioButton)
                    .addComponent(localjRadioButton)
                    .addComponent(affinejRadioButton))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        alignjPanelLayout.setVerticalGroup(
            alignjPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(alignjPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(globaljRadioButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(localjRadioButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(affinejRadioButton)
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        dnasqjPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("DNA Sequence"));

        dna1jLabel.setText("DNA sequence 1 :");

        dna2jLabel.setText("DNA sequence 2 :");

        dna1jTextField.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                dna1jTextFieldFocusLost(evt);
            }
        });

        dna2jTextField.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                dna2jTextFieldFocusLost(evt);
            }
        });

        javax.swing.GroupLayout dnasqjPanelLayout = new javax.swing.GroupLayout(dnasqjPanel);
        dnasqjPanel.setLayout(dnasqjPanelLayout);
        dnasqjPanelLayout.setHorizontalGroup(
            dnasqjPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(dnasqjPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(dnasqjPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addGroup(dnasqjPanelLayout.createSequentialGroup()
                        .addComponent(dna1jLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(dna1jTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 212, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(dnasqjPanelLayout.createSequentialGroup()
                        .addComponent(dna2jLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(dna2jTextField)))
                .addContainerGap(124, Short.MAX_VALUE))
        );
        dnasqjPanelLayout.setVerticalGroup(
            dnasqjPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(dnasqjPanelLayout.createSequentialGroup()
                .addGap(12, 12, 12)
                .addGroup(dnasqjPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(dna1jTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(dna1jLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(dnasqjPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(dna2jLabel)
                    .addComponent(dna2jTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        resultjPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Optimal Alignment"));

        resultjTextArea.setColumns(20);
        resultjTextArea.setRows(5);
        jScrollPane2.setViewportView(resultjTextArea);

        javax.swing.GroupLayout resultjPanelLayout = new javax.swing.GroupLayout(resultjPanel);
        resultjPanel.setLayout(resultjPanelLayout);
        resultjPanelLayout.setHorizontalGroup(
            resultjPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(resultjPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jScrollPane2, javax.swing.GroupLayout.DEFAULT_SIZE, 257, Short.MAX_VALUE)
                .addContainerGap())
        );
        resultjPanelLayout.setVerticalGroup(
            resultjPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(resultjPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jScrollPane2, javax.swing.GroupLayout.PREFERRED_SIZE, 115, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        penaltyjPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Penalty"));

        p1jLabel.setText("Gap penalty :");

        p2jLabel.setText("Affine gap penalty (introduction) :");

        p3jLabel.setText("Additional affine gap penalty :");

        javax.swing.GroupLayout penaltyjPanelLayout = new javax.swing.GroupLayout(penaltyjPanel);
        penaltyjPanel.setLayout(penaltyjPanelLayout);
        penaltyjPanelLayout.setHorizontalGroup(
            penaltyjPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(penaltyjPanelLayout.createSequentialGroup()
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addGroup(penaltyjPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(penaltyjPanelLayout.createSequentialGroup()
                        .addComponent(p1jLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(p1jTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 38, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(penaltyjPanelLayout.createSequentialGroup()
                        .addGroup(penaltyjPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addComponent(p3jLabel)
                            .addComponent(p2jLabel))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(penaltyjPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                            .addComponent(p2jTextField, javax.swing.GroupLayout.DEFAULT_SIZE, 38, Short.MAX_VALUE)
                            .addComponent(p3jTextField))))
                .addContainerGap())
        );
        penaltyjPanelLayout.setVerticalGroup(
            penaltyjPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(penaltyjPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(penaltyjPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(p1jLabel)
                    .addComponent(p1jTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(penaltyjPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(p2jLabel)
                    .addComponent(p2jTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(penaltyjPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(p3jLabel)
                    .addComponent(p3jTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(12, Short.MAX_VALUE))
        );

        scorejPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Score"));

        s1jTextField.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                s1jTextFieldFocusLost(evt);
            }
        });

        s2jTextField.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                s2jTextFieldFocusLost(evt);
            }
        });

        s3jTextField.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                s3jTextFieldFocusLost(evt);
            }
        });

        s1jLabel.setText("Score for a perfect match (i.e. (A,A), (C,C), (G,G), (T,T)) :");

        s2jLabel.setText("Score for a Purine-Purine or Pyrimidine-Pyrimidine substitution (i.e. (A,G), (G,A), (C,T), (T,C)) :");
        s2jLabel.setToolTipText("");

        s3jLabel.setText("Score for a Purine-Pyrimidine or Pyrimidine-Purine substitution :");

        javax.swing.GroupLayout scorejPanelLayout = new javax.swing.GroupLayout(scorejPanel);
        scorejPanel.setLayout(scorejPanelLayout);
        scorejPanelLayout.setHorizontalGroup(
            scorejPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(scorejPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(scorejPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(s2jLabel)
                    .addComponent(s1jLabel)
                    .addComponent(s3jLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(scorejPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(s1jTextField, javax.swing.GroupLayout.DEFAULT_SIZE, 50, Short.MAX_VALUE)
                    .addComponent(s2jTextField)
                    .addComponent(s3jTextField))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        scorejPanelLayout.setVerticalGroup(
            scorejPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(scorejPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(scorejPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(s1jTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(s1jLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(scorejPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(s2jLabel)
                    .addComponent(s2jTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(scorejPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(s3jLabel)
                    .addComponent(s3jTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        runjButton.setText("Run");
        runjButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                runjButtonActionPerformed(evt);
            }
        });

        clearjButton.setText("Clear");
        clearjButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                clearjButtonActionPerformed(evt);
            }
        });

        exitjButton.setText("Exit");
        exitjButton.setToolTipText("");
        exitjButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                exitjButtonActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout excutejPanelLayout = new javax.swing.GroupLayout(excutejPanel);
        excutejPanel.setLayout(excutejPanelLayout);
        excutejPanelLayout.setHorizontalGroup(
            excutejPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, excutejPanelLayout.createSequentialGroup()
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addComponent(clearjButton, javax.swing.GroupLayout.PREFERRED_SIZE, 88, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(runjButton, javax.swing.GroupLayout.PREFERRED_SIZE, 87, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(exitjButton, javax.swing.GroupLayout.PREFERRED_SIZE, 85, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(25, 25, 25))
        );
        excutejPanelLayout.setVerticalGroup(
            excutejPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, excutejPanelLayout.createSequentialGroup()
                .addGap(0, 9, Short.MAX_VALUE)
                .addGroup(excutejPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE, false)
                    .addComponent(clearjButton, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(runjButton, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(exitjButton))
                .addContainerGap())
        );

        tablejPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Score Table"));

        tablejTextArea.setColumns(20);
        tablejTextArea.setRows(5);
        jScrollPane3.setViewportView(tablejTextArea);

        javax.swing.GroupLayout tablejPanelLayout = new javax.swing.GroupLayout(tablejPanel);
        tablejPanel.setLayout(tablejPanelLayout);
        tablejPanelLayout.setHorizontalGroup(
            tablejPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(tablejPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jScrollPane3)
                .addContainerGap())
        );
        tablejPanelLayout.setVerticalGroup(
            tablejPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(tablejPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jScrollPane3, javax.swing.GroupLayout.PREFERRED_SIZE, 350, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                    .addComponent(scorejPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(dnasqjPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(tablejPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                            .addComponent(penaltyjPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(resultjPanel, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(alignjPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(excutejPanel, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.PREFERRED_SIZE, 0, Short.MAX_VALUE))))
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(scorejPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(penaltyjPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(alignjPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(resultjPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(12, 12, 12)
                        .addComponent(excutejPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(dnasqjPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(tablejPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
                .addContainerGap())
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    
    
    private void exitjButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_exitjButtonActionPerformed
        System.exit(0);
    }//GEN-LAST:event_exitjButtonActionPerformed

    private void clearjButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_clearjButtonActionPerformed
        p1jTextField.setText("");
        p2jTextField.setText("");
        p3jTextField.setText("");
        s1jTextField.setText("");
        s2jTextField.setText("");
        s3jTextField.setText("");
        dna1jTextField.setText("");
        dna2jTextField.setText("");
        resultjTextArea.setText("");
    }//GEN-LAST:event_clearjButtonActionPerformed

    private void runjButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_runjButtonActionPerformed
        int s1 = Integer.parseInt(s1jTextField.getText());
        int s2 = Integer.parseInt(s2jTextField.getText());
        int s3 = Integer.parseInt(s3jTextField.getText());
                
        int p1 = Integer.parseInt(p1jTextField.getText());
        int p2 = Integer.parseInt(p2jTextField.getText());
        int p3 = Integer.parseInt(p3jTextField.getText());
        
        String sq1 = dna1jTextField.getText();
        String sq2 = dna2jTextField.getText();
        
        /*-----------------------------------------------------------------*/
        // Error Message when input is invalid
        try {
            Integer.parseInt(s1jTextField.getText());
            Integer.parseInt(s2jTextField.getText());
            Integer.parseInt(s3jTextField.getText());
            Integer.parseInt(p1jTextField.getText());
            Integer.parseInt(p2jTextField.getText());
            Integer.parseInt(p3jTextField.getText());
        }
        catch (Exception e)
        {
            JOptionPane.showMessageDialog(null, "The input is invalid", "Error", JOptionPane.ERROR_MESSAGE);
            return;
        }
        
        if(Integer.parseInt(s1jTextField.getText()) < 0) 
        {
                JOptionPane.showMessageDialog(null, "The input is invalid", "Error", JOptionPane.ERROR_MESSAGE);  
        }
        
        if (Integer.parseInt(s1jTextField.getText()) < Integer.parseInt(s2jTextField.getText()))
        {
            JOptionPane.showMessageDialog(null, "The input is invalid", "Error", JOptionPane.ERROR_MESSAGE);
        }
        if (Integer.parseInt(s1jTextField.getText()) < Integer.parseInt(s2jTextField.getText()))
        {
            JOptionPane.showMessageDialog(null, "The input is invalid", "Error", JOptionPane.ERROR_MESSAGE);
        }
        
        if (Integer.parseInt(s2jTextField.getText()) < Integer.parseInt(s3jTextField.getText()))
        {
            JOptionPane.showMessageDialog(null, "The input is invalid", "Error", JOptionPane.ERROR_MESSAGE);
        }
        if (Integer.parseInt(s2jTextField.getText()) < Integer.parseInt(s3jTextField.getText()))
        {
            JOptionPane.showMessageDialog(null, "The input is invalid", "Error", JOptionPane.ERROR_MESSAGE);
        }
        
        if (!sq1.matches("[ACGT]*"))
        {
            JOptionPane.showMessageDialog(null, "Invalid DNA sequence", "Error", JOptionPane.ERROR_MESSAGE);
        }
        if (!sq2.matches("[ACGT]*"))
        {
            JOptionPane.showMessageDialog(null, "Invalid DNA sequence", "Error", JOptionPane.ERROR_MESSAGE);
        }
        
        
        /*-----------------------------------------------------------------*/
        // Methods for alignment
        
        if (global)
        {
            String globalResult = method.global(s1, s2, s3, p1, sq1, sq2);
            resultjTextArea.setText(globalResult);
            
            tablejTextArea.setText(method.g.toString());
        }
        else if(local)
        {
            String localResult = method.local(s1, s2, s3, p1, sq1, sq2);
            resultjTextArea.setText(localResult);
            
            tablejTextArea.setText(method.l.toString());
        }
        else if(affine)
        {
            String affineResult = method.affine(s1, s2, s3, p2, p3, sq1, sq2);
            resultjTextArea.setText(affineResult);
            
            tablejTextArea.setText(method.af.toString());
        } 
    }//GEN-LAST:event_runjButtonActionPerformed

    private void globaljRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_globaljRadioButtonActionPerformed
        this.global = true;
        this.local = false;
        this.affine = false;
        localjRadioButton.setSelected(false);
        affinejRadioButton.setSelected(false);
    }//GEN-LAST:event_globaljRadioButtonActionPerformed

    private void localjRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_localjRadioButtonActionPerformed
        this.local = true;
        this.global = false;
        this.affine = false;
        globaljRadioButton.setSelected(false);
        affinejRadioButton.setSelected(false);
    }//GEN-LAST:event_localjRadioButtonActionPerformed

    private void affinejRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_affinejRadioButtonActionPerformed
        this.affine = true;
        this.local = false;
        this.global = false;
        localjRadioButton.setSelected(false);
        globaljRadioButton.setSelected(false);
    }//GEN-LAST:event_affinejRadioButtonActionPerformed

    private void s1jTextFieldFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_s1jTextFieldFocusLost
        try {
            Integer.parseInt(s1jTextField.getText());
        }
        catch (Exception e)
        {
            JOptionPane.showMessageDialog(null, "Invalid Input", "Error", JOptionPane.ERROR_MESSAGE);
            return;
        }
        
        if(Integer.parseInt(s1jTextField.getText()) < 0) 
        {
                JOptionPane.showMessageDialog(null, "The score must be positive", "Error", JOptionPane.ERROR_MESSAGE);  
        }
        
        if (Integer.parseInt(s1jTextField.getText()) < Integer.parseInt(s2jTextField.getText()))
        {
            JOptionPane.showMessageDialog(null, "The score must be larger than the score for the same bases", "Error", JOptionPane.ERROR_MESSAGE);
        }
    }//GEN-LAST:event_s1jTextFieldFocusLost

    private void s2jTextFieldFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_s2jTextFieldFocusLost
        try {
            Integer.parseInt(s2jTextField.getText());
        }
        catch (Exception e)
        {
            JOptionPane.showMessageDialog(null, "Invalid Input", "Error", JOptionPane.ERROR_MESSAGE);
            return;
        }
        
        if (Integer.parseInt(s1jTextField.getText()) < Integer.parseInt(s2jTextField.getText()))
        {
            JOptionPane.showMessageDialog(null, "The score must be less than the score for perfect match", "Error", JOptionPane.ERROR_MESSAGE);
        }
        
        if (Integer.parseInt(s2jTextField.getText()) < Integer.parseInt(s3jTextField.getText()))
        {
            JOptionPane.showMessageDialog(null, "The score must be larger than the score for the different bases", "Error", JOptionPane.ERROR_MESSAGE);
        }
    }//GEN-LAST:event_s2jTextFieldFocusLost

    private void s3jTextFieldFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_s3jTextFieldFocusLost
        try {
            Integer.parseInt(s3jTextField.getText());
        }
        catch (Exception e)
        {
            JOptionPane.showMessageDialog(null, "Invalid Input", "Error", JOptionPane.ERROR_MESSAGE);
            return;
        }
        
        if (Integer.parseInt(s2jTextField.getText()) < Integer.parseInt(s3jTextField.getText()))
        {
            JOptionPane.showMessageDialog(null, "The score must be less than the score for the same bases", "Error", JOptionPane.ERROR_MESSAGE);
        }
    }//GEN-LAST:event_s3jTextFieldFocusLost

    private void dna1jTextFieldFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_dna1jTextFieldFocusLost
        String sq1 = dna1jTextField.getText();
        if (!sq1.matches("[ACGT]*"))
        {
            JOptionPane.showMessageDialog(null, "Invalid DNA sequence", "Error", JOptionPane.ERROR_MESSAGE);
        }
    }//GEN-LAST:event_dna1jTextFieldFocusLost

    private void dna2jTextFieldFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_dna2jTextFieldFocusLost
         String sq2 = dna2jTextField.getText();
        if (!sq2.matches("[ACGT]*"))
        {
            JOptionPane.showMessageDialog(null, "Invalid DNA sequence", "Error", JOptionPane.ERROR_MESSAGE);
        }
    }//GEN-LAST:event_dna2jTextFieldFocusLost

    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        /* Set the Nimbus look and feel */
        //<editor-fold defaultstate="collapsed" desc=" Look and feel setting code (optional) ">
        /* If Nimbus (introduced in Java SE 6) is not available, stay with the default look and feel.
         * For details see http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html 
         */
        try {
            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
                if ("Nimbus".equals(info.getName())) {
                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
        } catch (ClassNotFoundException ex) {
            java.util.logging.Logger.getLogger(SAUI.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (InstantiationException ex) {
            java.util.logging.Logger.getLogger(SAUI.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (IllegalAccessException ex) {
            java.util.logging.Logger.getLogger(SAUI.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(SAUI.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>

        /* Create and display the form */
        java.awt.EventQueue.invokeLater(new Runnable() {
            public void run() {
                new SAUI().setVisible(true);
            }
        });
    }

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JRadioButton affinejRadioButton;
    private javax.swing.JPanel alignjPanel;
    private javax.swing.JButton clearjButton;
    private javax.swing.JLabel dna1jLabel;
    private javax.swing.JTextField dna1jTextField;
    private javax.swing.JLabel dna2jLabel;
    private javax.swing.JTextField dna2jTextField;
    private javax.swing.JPanel dnasqjPanel;
    private javax.swing.JPanel excutejPanel;
    private javax.swing.JButton exitjButton;
    private javax.swing.JRadioButton globaljRadioButton;
    private javax.swing.JScrollPane jScrollPane2;
    private javax.swing.JScrollPane jScrollPane3;
    private javax.swing.JRadioButton localjRadioButton;
    private javax.swing.JLabel p1jLabel;
    private javax.swing.JTextField p1jTextField;
    private javax.swing.JLabel p2jLabel;
    private javax.swing.JTextField p2jTextField;
    private javax.swing.JLabel p3jLabel;
    private javax.swing.JTextField p3jTextField;
    private javax.swing.JPanel penaltyjPanel;
    private javax.swing.JPanel resultjPanel;
    private javax.swing.JTextArea resultjTextArea;
    private javax.swing.JButton runjButton;
    private javax.swing.JLabel s1jLabel;
    private javax.swing.JTextField s1jTextField;
    private javax.swing.JLabel s2jLabel;
    private javax.swing.JTextField s2jTextField;
    private javax.swing.JLabel s3jLabel;
    private javax.swing.JTextField s3jTextField;
    private javax.swing.JPanel scorejPanel;
    private javax.swing.JPanel tablejPanel;
    private javax.swing.JTextArea tablejTextArea;
    // End of variables declaration//GEN-END:variables

}