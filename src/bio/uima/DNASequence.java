

/* First created by JCasGen Sat Mar 30 12:47:00 EDT 2013 */
package bio.uima;

import org.apache.uima.jcas.JCas; 
import org.apache.uima.jcas.JCasRegistry;
import org.apache.uima.jcas.cas.TOP_Type;

import org.apache.uima.jcas.tcas.Annotation;


/** 
 * Updated by JCasGen Thu Apr 04 23:16:27 EDT 2013
 * XML source: /home/farhang/workspace/BioUIMA/desc/BioDescriptor.xml
 * @generated */
public class DNASequence extends Annotation {
  /** @generated
   * @ordered 
   */
  @SuppressWarnings ("hiding")
  public final static int typeIndexID = JCasRegistry.register(DNASequence.class);
  /** @generated
   * @ordered 
   */
  @SuppressWarnings ("hiding")
  public final static int type = typeIndexID;
  /** @generated  */
  @Override
  public              int getTypeIndexID() {return typeIndexID;}
 
  /** Never called.  Disable default constructor
   * @generated */
  protected DNASequence() {/* intentionally empty block */}
    
  /** Internal - constructor used by generator 
   * @generated */
  public DNASequence(int addr, TOP_Type type) {
    super(addr, type);
    readObject();
  }
  
  /** @generated */
  public DNASequence(JCas jcas) {
    super(jcas);
    readObject();   
  } 

  /** @generated */  
  public DNASequence(JCas jcas, int begin, int end) {
    super(jcas);
    setBegin(begin);
    setEnd(end);
    readObject();
  }   

  /** <!-- begin-user-doc -->
    * Write your own initialization here
    * <!-- end-user-doc -->
  @generated modifiable */
  private void readObject() {/*default - does nothing empty block */}
     
 
    
  //*--------------*
  //* Feature: value

  /** getter for value - gets Stores the DNA Sequence string
   * @generated */
  public String getValue() {
    if (DNASequence_Type.featOkTst && ((DNASequence_Type)jcasType).casFeat_value == null)
      jcasType.jcas.throwFeatMissing("value", "bio.uima.DNASequence");
    return jcasType.ll_cas.ll_getStringValue(addr, ((DNASequence_Type)jcasType).casFeatCode_value);}
    
  /** setter for value - sets Stores the DNA Sequence string 
   * @generated */
  public void setValue(String v) {
    if (DNASequence_Type.featOkTst && ((DNASequence_Type)jcasType).casFeat_value == null)
      jcasType.jcas.throwFeatMissing("value", "bio.uima.DNASequence");
    jcasType.ll_cas.ll_setStringValue(addr, ((DNASequence_Type)jcasType).casFeatCode_value, v);}    
  }

    