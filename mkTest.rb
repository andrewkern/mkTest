# mkTest.rb
#
# Andrew Kern
#
##############
## Ruby implementation of the McDonald-Kreitman test
#  and a whole lot more. 
#
#
#
#
#
# Supplied in this script is my complete Ruby popgen toolkit, so a totally bloated script but easy to give to friends and colleagues to use
# feel free to use the code below, alter it, and distribute it how you see fit but please give me some credit somewhere along the way!

include Math


module ArrayADK

	def checkArray		#returns edited self or nil 
		temp = self.delete_if{ | x | x == "NA" }
		if temp.empty?
			return nil
		else
			return temp
		end
	end
	def inject(n)
		each do | value |
			n = yield(n, value)
		end
		n
	end
	def sum(initial = 0)
		inject(initial){ | n, value | n + value }
	end
	def max
		temp = self.checkArray
		if temp
			max = 0.0
			temp.each{ | x | 
					if x > max
						max = x
					end
					}
			return max
		else
			return "NA"
		end
	end
	def product(initial = 1)
		inject(initial){ | n, value | n * value }
	end
	def sampleMean
		temp = self.delete_if{ | x | x.class == String }
		if temp.empty?
			return "NA"
		else
			return (temp.sum).to_f/(temp.size).to_f
		end
	end
        def weightedMean(weightArray)
            numSum = denSum = 0.0
            0.upto(self.size - 1){ | i |
                numSum += self[i].to_f * weightArray[i].to_f
                denSum += weightArray[i].to_f
            }
            return(numSum / denSum)
        end
                
	def sampleVariance
		temp = self.delete_if{ | x | x.class == String }
		mean = temp.sampleMean
		count = 0
		temp.each{ | x |
				count += (x - mean) * (x - mean)
				}
		count/(temp.size - 1)
	end
	def stdev
		temp = self.delete_if{ | x | x.class == String }
		mean = temp.sampleMean
		count = 0
		temp.each{ | x |
				count += (x - mean) * (x - mean)
				}
		return Math.sqrt(count/temp.size)
	end
	
	def occurrencesOf(anObject)
		count = 0
		self.each{ | x | 
			if x == anObject
				count += 1
			end }
		return count
	end
	def pearsonCorrelation(secondArray)
		if self.size == secondArray.size
			naArray = Array.new
			if self.include?("NA") or secondArray.include?("NA")
				self.each_index{ | i | if self[i] == "NA" or secondArray[i] == "NA"
											naArray<<i
										end
									}
			end
			naArray.reverse.each{ | x | self.delete_at(x)
										secondArray.delete_at(x)
									}
			
			
			temp1 = self.checkArray
			temp2 = secondArray.checkArray
			
			mean1 = temp1.sampleMean
			mean2 = temp2.sampleMean
			
			stdev1 = temp1.stdev
			stdev2 = temp2.stdev
			sum = 0
			
			0.upto(temp1.size - 1){ | i |
				sum += ((temp1[i] - mean1) / stdev1) * ((temp2[i] - mean2) / stdev2)
				}
			return sum / temp1.size.to_f
			
			
			
		else
			print "Error: Arrays different lengths\n"
			exit()
		end
	end
	
	
	def distancesBetween(anObject)   #returns an array of index distances between anObject in current array
		dists = Array.new
		positions = Array.new
		self.each_index{ | i |
				if (self[i] == anObject)
					positions.push(i)
				end
				}
		positions.reverse!
		positions.each_index { | i |
					if i < (positions.size - 1)
						dists.push( positions[i] - positions[i + 1] + 1)
					end
					}
		return dists
	end
	
	def distancesBetweenTwo(anObject, anotherObject)   #returns an array of index distances between anObject in current array
		dists = Array.new
		positions = Array.new
		current = anObject
		self.each_index{ | i |
				if (self[i].to_s == current)
					positions.push(i)
					if current == anObject
						current = anotherObject
					else
						current = anObject
					end
				end
				}
		if positions.empty?
			print "distancesBetween(anObject): no match in to object!\n"
			exit
		end
		positions.reverse!
		positions.each_index { | i |
					if i < (positions.size - 1)
						dists.push( positions[i] - positions[i + 1])
					end
					}
		return dists
	end
	
	def randomArray
		rSum = self[0] + self[1]  #returns a random array of length 2
		cSum = self[0] + self[2]
		min = [rSum, cSum].min
		new = Array.new(2,0)
		new[0] = rand(min)
		new[1] = rSum - new[0] 
		return new
	end
	
	def randomTable   #takes an array with 4 entries, returns random array conditional on marginals
		new = Array.new
		col1 = self[0] + self[2]
		col2 = self[1] + self[3]
		randRow = self.randomArray
		new.push(randRow[0])
		new.push(randRow[1])
		new.push(col1 - randRow[0])
		new.push(col2 - randRow[1])
		return new
	end
	
	def monteCarloContingency		#take an array with 4 entries
		chiArray = Array.new
		10000.times { 
			chiArray.push(self.randomTable.chisquare)
			}
		chiArray.sort!
		
		testStat = self.chisquare
		chiArray.each_index{ | i |
					if chiArray[i] >= testStat
						return (1 - (i.to_f/10000.0))
					end
					}
	
	end
	
end




# CodingSequence Module - This gives methods to the SequenceMatrix class which describe coding sequence
##
##

module CodingSequence

  #coding sequences
  def codingSites
    cr = Array.new
    self.codingRegions.each{ | x | cr.push(x - 1) }
    oc = Array.new
    exonNumber = cr.size / 2
    exonNumber.times{ | i |
      b = cr.shift
      e = cr.shift
      b.upto(e){ | site |
        oc.push(site) }
      }
      return oc	
    end

    def intronSites			# this assumes only those sites surrounded by coding sequence are introns
      is = Array.new
      cs = self.codingSites
      if cs.empty?
        0.upto(self.matrix[1].size - 1){ | i | is.push(i) }
      end
      (self.codingRegions.first - 1).upto(self.codingRegions.last - 1){ | i | is.push(i) }
      cs.each{ | x | is.delete(x)}
      return is
    end

    def setCodons
      trips = Array.new
      oc = Array.new
      pad = self.readingFrame - 1

      self.matrix.each{ | anAllele |
        temp = String.new("")
        self.codingSites.each{ | c |
          temp << anAllele[c,1].to_s 
        }
        oc.push(temp)}
        oc.each{ | cds |
          pad.times{
            cds.insert(0,"-") }
          }
          codonNumber = (oc[0].size / 3.0 ).floor
          top = oc[0].size
          lastCompleteCodonBase = codonNumber * 3
          if (top - lastCompleteCodonBase) != 0
            oc.each{ | cds |
              (3 - (top - lastCompleteCodonBase)).times{ cds << ("-") }
            }
          end

          oc.each{ | cds |
            temp = cds.scan(/.../)
            trips.push(temp)
          }
          self.codons = trips	
          return self
        end
        def codonSet(aCodonIndex)
          cSet = Array.new
          self.codons.each{ | each | cSet.push(each[aCodonIndex])}
          return cSet.uniq
        end

        def codonSetClean(aCodonIndex)
          cSet = codonSet(aCodonIndex)
          cSet.reject!{ | x | x.include?("N") or x.include?("-") }
          return cSet
        end

        def aminoSet(aCodonIndex)
          aSet = Array.new
          self.codons.each{ | each | aSet.push(self.geneticCode[each[aCodonIndex]]) }
          return aSet.uniq
        end

        def puSet(aCodonIndex)
          c = self.codonSet(aCodonIndex)
          aSet = Array.new
          edit = c.reject{ | each | each.include?("N") or each.include?("-") }
          edit.each{ | each | aSet.push(self.puDict[each]) }
          return aSet.uniq
        end

        def aminoCodonSet(aSet)
          edit = aSet.reject{ | each | each.include?("N") or each.include?("-") }
          anArray = Array.new
          edit.each{ | each | 
            anArray << self.geneticCode[each.upcase]
          }
          return anArray.uniq
        end

        def codonSegSites(aCodonIndex)
          cSet = Array.new
          self.codons.each{ | anArray |
            cSet.push(anArray[aCodonIndex])
          }
          cSet.uniq!
          setHold = Array.new(3,Array.new)

          sum = 0
          0.upto(2){ | i | 
            test = Array.new
            cSet.each{ | codon | test.push(codon[i,1])}
            test.uniq!
            if test.size > 1
              sum += 1
            end }

            return sum
          end

          def codonSegSitesSet(aCodonSet)						#returns number of sites segregating in a codon set
            sum = 0
            0.upto(2){ | c |
              temp = Array.new
              aCodonSet.each{ | eachCodon |
                temp << eachCodon[c,1]
              }
              if temp.include?("N") or temp.include?("-")
                return 0
              else
                if temp.uniq.size > 1
                  sum += 1
                end
              end
            }
            return sum
          end

          def codonSegSitesPositions(aCodonSet)				# returns an Array of positions
            pos = Array.new
            sum = 0
            0.upto(2){ | c |
              temp = Array.new
              aCodonSet.each{ | eachCodon |
                temp << eachCodon[c,1]
              }
              if temp.include?("N") or temp.include?("-")
                return Array.new
              else
                if temp.uniq.size > 1
                  pos << c
                end
              end
            }
            return pos
          end

          def codonSegSitesPositionDict  						#returns array of segSites in both codon and nuc. indexes
            pos = Hash.new
            segSites = self.segSitesLocationsAllSites
            codingSites = self.codingSites
            segs = segSites.find_all{ | x | codingSites.include?(x) }
            if segs.empty?
              return nil
            else
              codonPos = Array.new
              segs.each{ | each | codonPos << (((codingSites.index(each)).to_f) / 3) }
              pos["codonSegSites"] = codonPos
              pos["locusSegSites"] = segs
              return pos
            end
          end

          def codonIndexToCodingSequenceIndex(anIndex)  						#returns first nucleotide index of codon
            codingSites = self.codingSites
            nucIndex = anIndex * 3
            return codingSites[nucIndex]
          end

          def shortestPath2Codons(aCodonSet)					# takes a set of codons and returns an array indicating a path

            aminoSet = self.aminoCodonSet(aCodonSet)
            pos = self.codonSegSitesPositions(aCodonSet)
            segSites = pos.size
            paths = Array.new

            if segSites == 1
              path = Array.new
              if aminoSet.size > 1
                path << ["R",pos[0].to_i]
              else
                path << ["S",pos[0].to_i]
              end
              paths << path
            else
              if segSites == 2
                perms = [[pos[0],pos[1]],[pos[1],pos[0]]]    #permutations for swaps
                paths = Array.new
                tempString = String.new("")
                tempString << aCodonSet[0]
                perms.each{ | eachArray |

                  path = Array.new
                  tempString = String.new("")
                  tempString << aCodonSet[0]
                  eachArray.each{ | i |
                    tempAmino1 = String.new("")
                    tempAmino1 << self.geneticCode[tempString]
                    tempString[i,1] = aCodonSet[1][i,1]
                    tempAmino2 = String.new("")
                    tempAmino2 << self.geneticCode[tempString]
                    if tempAmino2 == "*"
                      path << ["*",i]
                    end	
                    if tempAmino1 == tempAmino2
                      path << ["S", i]
                    else
                      path << ["R", i]
                    end		
                  }
                  paths << path }
                else	#3 segSites!!!
                  perms = [[0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,1,0],[2,0,1]]    #permutations for swaps
                  paths = Array.new
                  perms.each{ | eachArray |
                    path = Array.new
                    tempString = String.new("")
                    tempString << aCodonSet[0]
                    eachArray.each{ | i |
                      tempAmino1 = String.new("")
                      tempAmino1 << self.geneticCode[tempString]
                      tempString[i,1] = aCodonSet[1][i,1]
                      tempAmino2 = String.new("")
                      tempAmino2 << self.geneticCode[tempString]
                      if tempAmino2 == "*"
                        path << ["*",i]
                      end	
                      if tempAmino1 == tempAmino2
                        path << ["S", i]
                      else
                        path << ["R", i]
                      end		
                    }
                    paths << path }
                  end
                end		
                #set up the path length variable
                short = 100
                shortPath = nil
                noStopsPaths = paths.delete_if{ | x | x.first.first == "*" }
                #go through paths, find shortest with respect to Replacements
                noStopsPaths.each{ | aPath |
                  rLength = 0
                  aPath.each{ | eachArray | if eachArray.include?("R")
                    rLength += 1 
                  end
                }
                if short > rLength
                  short = rLength
                  shortPath = aPath
                else
                  if shortPath.size > aPath.size		# more replacements but less changes overall
                    shortPath = aPath
                  end
                end
              }
              return shortPath

            end
            def lookupPath2Codons(aCodonSet)
              if aCodonSet.include?("TAA") or aCodonSet.include?("TAG") or aCodonSet.include?("TGA")
                return nil
              end
              string = self.codonDists[aCodonSet.first][aCodonSet.last]
              path = Array.new
              steps = string.length / 2
              0.upto(steps - 1){ | i | path << [string[(2*i),1],string[(2*i)+1,1]] }
              return path
            end

            def shortestPath3Codons(aCodonSet)
              segSites = self.codonSegSitesSet(aCodonSet)
              aminoSet = self.aminoCodonSet(aCodonSet)
              pos = self.codonSegSitesPositions(aCodonSet)
              paths = Array.new
              shortPath = nil

              if segSites == 1
                if aminoSet.size > 1
                  path = Array.new
                  r = aminoSet.size - 1
                  s = aCodonSet.size - 1 - r
                  r.times{ path << ["R",pos[0]] }
                  s.times{ path << ["S",pos[0]] }
                  paths << path
                else
                  path = Array.new
                  s = aCodonSet.size - 1
                  s.times{ path << ["S",pos[0]] }
                  paths << path
                end
              else
                if segSites == 2
                  if aCodonSet.size == 3
                    #set up all possible paths, use 2 codon paths

                    codonPaths = Array.new
                    codonPaths << [[aCodonSet[0],aCodonSet[1]],[aCodonSet[1],aCodonSet[2]]]
                    codonPaths << [[aCodonSet[0],aCodonSet[2]],[aCodonSet[1],aCodonSet[2]]]
                    codonPaths << [[aCodonSet[0],aCodonSet[2]],[aCodonSet[0],aCodonSet[1]]]
                    codonPaths.each{ | eachPath |
                      temp = Array.new
                      eachPath.each{ | testArray |
                        if self.lookupPath2Codons(testArray).nil?
                          print testArray,"\n"
                        end
                        self.lookupPath2Codons(testArray).each{ | x | temp << x }
                      }
                      paths << temp
                    }
                  else
                    if aCodonSet.size == 4
                      #set up all possible paths, use 2 codon paths
                      codonPaths = Array.new
                      codonPaths << [[aCodonSet[0],aCodonSet[1]],[aCodonSet[1],aCodonSet[2]],[aCodonSet[2],aCodonSet[3]]]				
                      codonPaths << [[aCodonSet[0],aCodonSet[1]],[aCodonSet[1],aCodonSet[3]],[aCodonSet[3],aCodonSet[2]]]
                      codonPaths << [[aCodonSet[0],aCodonSet[2]],[aCodonSet[2],aCodonSet[1]],[aCodonSet[1],aCodonSet[3]]]
                      codonPaths << [[aCodonSet[0],aCodonSet[2]],[aCodonSet[2],aCodonSet[3]],[aCodonSet[3],aCodonSet[1]]]
                      codonPaths << [[aCodonSet[0],aCodonSet[3]],[aCodonSet[3],aCodonSet[2]],[aCodonSet[2],aCodonSet[1]]]
                      codonPaths << [[aCodonSet[0],aCodonSet[3]],[aCodonSet[3],aCodonSet[1]],[aCodonSet[1],aCodonSet[2]]]
                      codonPaths << [[aCodonSet[1],aCodonSet[0]],[aCodonSet[0],aCodonSet[2]],[aCodonSet[2],aCodonSet[3]]]
                      codonPaths << [[aCodonSet[1],aCodonSet[3]],[aCodonSet[3],aCodonSet[0]],[aCodonSet[0],aCodonSet[2]]]
                      codonPaths << [[aCodonSet[1],aCodonSet[3]],[aCodonSet[1],aCodonSet[0]],[aCodonSet[0],aCodonSet[2]]]
                      codonPaths << [[aCodonSet[1],aCodonSet[2]],[aCodonSet[2],aCodonSet[0]],[aCodonSet[0],aCodonSet[3]]]
                      codonPaths << [[aCodonSet[1],aCodonSet[0]],[aCodonSet[0],aCodonSet[3]],[aCodonSet[3],aCodonSet[2]]]
                      codonPaths << [[aCodonSet[3],aCodonSet[0]],[aCodonSet[0],aCodonSet[1]],[aCodonSet[1],aCodonSet[2]]]
                      codonPaths << [[aCodonSet[0],aCodonSet[1]],[aCodonSet[0],aCodonSet[3]],[aCodonSet[0],aCodonSet[2]]]
                      codonPaths << [[aCodonSet[1],aCodonSet[0]],[aCodonSet[1],aCodonSet[2]],[aCodonSet[1],aCodonSet[3]]]
                      codonPaths << [[aCodonSet[3],aCodonSet[2]],[aCodonSet[3],aCodonSet[1]],[aCodonSet[3],aCodonSet[0]]]
                      codonPaths << [[aCodonSet[2],aCodonSet[0]],[aCodonSet[2],aCodonSet[1]],[aCodonSet[2],aCodonSet[3]]]
                      codonPaths.each{ | eachPath |
                        temp = Array.new
                        eachPath.each{ | testArray |
                          self.lookupPath2Codons(testArray).each{ | x | temp << x }
                        }
                        paths << temp
                      }
                    else
                      return ["C",pos[0]]
                    end
                  end
                else
                  return ["C",pos[0]]
                end
              end	

              #set up the path length variable
              short = 100
              shortPath = nil
              paths = paths.reject{ | x | x.include?("*") }

              #go through paths, find shortest with respect to Replacements
              paths.each{ | aPath |
                rLength = 0
                aPath.each{ | eachArray | if eachArray.include?("R")
                  rLength += 1 
                end
              }
              if short > rLength
                short = rLength
                shortPath = aPath
              else
                if shortPath.size > aPath.size		# more replacements but less changes overall
                  shortPath = aPath
                end
              end
            }
            return shortPath
          end

          def allPaths2Codons(aCodonSet)		# takes a set of codons and returns an array of paths
            #initialize arrays
            aminoSet = self.aminoCodonSet(aCodonSet)
            pos = self.codonSegSitesPositions(aCodonSet)
            segSites = pos.size
            paths = Array.new

            #1 segSite?
            if segSites == 1
              path = Array.new
              if aminoSet.size > 1  #is it a aa replacement?
                path << ["R",pos[0].to_i]
              else
                path << ["S",pos[0].to_i]
              end
              paths << path
            else
              #2 segsites?
              if segSites == 2
                perms = [[pos[0],pos[1]],[pos[1],pos[0]]]    #permutations for swaps
                tempString = String.new("")
                tempString << aCodonSet[0]
                perms.each{ | eachArray |
                  path = Array.new
                  tempString = String.new("")
                  tempString << aCodonSet[0]
                  eachArray.each{ | i |
                    tempAmino1 = String.new("")
                    tempAmino1 << self.geneticCode[tempString]
                    tempString[i,1] = aCodonSet[1][i,1]
                    tempAmino2 = String.new("")
                    tempAmino2 << self.geneticCode[tempString]
                    if tempAmino2 == "*"
                      path << ["*",i]
                    end	
                    if tempAmino1 == tempAmino2
                      path << ["S", i]
                    else
                      path << ["R", i]
                    end		
                  }
                  paths << path }
                else	#3 segSites!!!
                  perms = [[0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,1,0],[2,0,1]]    #permutations for swaps
                  paths = Array.new
                  perms.each{ | eachArray |
                    path = Array.new
                    tempString = String.new("")
                    tempString << aCodonSet[0]
                    eachArray.each{ | i |
                      tempAmino1 = String.new("")
                      tempAmino1 << self.geneticCode[tempString]
                      tempString[i,1] = aCodonSet[1][i,1]
                      tempAmino2 = String.new("")
                      tempAmino2 << self.geneticCode[tempString]
                      if tempAmino2 == "*"
                        path << ["*",i]
                      end	
                      if tempAmino1 == tempAmino2
                        path << ["S", i]
                      else
                        path << ["R", i]
                      end		
                    }
                    paths << path 
                  }
                end
              end
              return paths
            end

            def silReplMutations		#returns a dictionary
              subs = Hash.new
              last = self.codons[0].length - 1
              ["replacements","silents","missingData","complex"].each{ | x | subs[x] = Array.new}
              codingSites = self.codingSites
              dict = self.codonSegSitesPositionDict
              if dict
                segLocus = dict["locusSegSites"]
                segCodons = dict["codonSegSites"]
                #reduce segCodons
                segCodons.collect!{| x | x.floor }.uniq!
                #go through reduced set and do stuff
                segCodons.each_with_index{ | site, index |
                  if site.floor != last
                    #  print site,"\n"
                    set = self.codonSet(site)
                    mis = 0

                    #is there ambiguous data? if yes cleanup

                    set.each{ | each | if each.include?("N") or each.include?("-")
                      mis += 1
                    end
                  }
                  if mis > 0
                    subs["missingData"] << segLocus[index]
                    set = set.reject{ | x | x.include?("N") or x.include?("-") }
                  end

                  #how many states?

                  if set.uniq.size == 2				
                    shortPath = lookupPath2Codons(set.uniq)
                    if shortPath
                      shortPath.each{ | each |
                        if each.include?("R")
                          subs["replacements"] << (site * 3) + each[1].to_i #segLocus[index]
                        else
                          subs["silents"] <<  (site * 3) + each[1].to_i # segLocus[index]
                        end
                      }
                    end
                  else
                    shortPath = self.shortestPath3Codons(set.uniq)
                    if shortPath
                      if shortPath.include?("C")
                        subs["complex"] << (site * 3) #segLocus[index]
                      else
                        shortPath.each{ | each | 
                          if each.include?("R")
                            subs["replacements"] << (site * 3) + each[1].to_i #segLocus[index]
                          else
                            subs["silents"] << (site * 3) + each[1].to_i #segLocus[index]
                          end
                        }
                      end
                    end
                  end
                end
              }
            end
            return subs
          end

          def silReplDifsAllPaths2Seqs(seqIndex1, seqIndex2)  #this returns an array (nonSyn/nsSites, Syn/sSites) as estimated by Nei & Gojobori alg. 1

            codons1 = self.codons[seqIndex1]
            codons2 = self.codons[seqIndex2]

            #initialize arrays
            nsCount = 0
            sCount = 0
            nsSites = 0
            sSites = 0
            #go through codons, count up number of diffs in each codon and their positions
            codons1.each_index{ | i |
              set = Array.new 
              set << codons1[i]
              set << codons2[i]
              #get rid of any ambiguous codons
              set.reject!{ | x | x.include?("N") or x.include?("X") or self.geneticCode[x] == "*" }
              if set.size > 1
                #tally up sites along the way
                nsSites += self.silentSiteDict[set[0]][1]
                nsSites += self.silentSiteDict[set[1]][1]
                sSites += self.silentSiteDict[set[0]][0]
                sSites += self.silentSiteDict[set[1]][0]
              end
              #more than one codon state?
              if set.uniq.size > 1
                #get allPaths between codons
                paths = self.allPaths2Codons(set)
                tempRCount = 0
                tempSCount = 0
                #clean out paths which go through stop codons
                paths.reject!{ | x | x[0] == "*"}
                #go through paths, count stuff
                paths.each{ | aPath |
                  #go through path steps, count em
                  aPath.each{ | aStep |
                    if aStep[0] == "R"
                      tempRCount += 1
                    else
                      tempSCount += 1
                    end
                  }
                }
                nsCount += tempRCount.to_f / paths.size
                sCount += tempSCount.to_f/ paths.size
              end
            }
            if (nsSites == 0 and sSites == 0)
              return [nil,nil]
            else
              return [nsCount / (nsSites.to_f / 2), sCount / (sSites.to_f / 2)]    
            end
          end


          #accessing info

          def setGeneticCode(aCodeTable)

            if aCodeTable == "standard"
              gc = { "GCT" => "A",
                "GCC" => "A",
                "GCA" => "A",
                "GCG" => "A",
                "TGT" => "C",
                "TGC" => "C",
                "GAT" => "D",
                "GAC" => "D",
                "GAA" => "E",
                "GAG" => "E",
                "TTT" => "F",
                "TTC" => "F",
                "GGT" => "G",
                "GGC" => "G",
                "GGA" => "G",
                "GGG" => "G",
                "CAT" => "H",
                "CAC" => "H",
                "ATT" => "I",
                "ATC" => "I",
                "ATA" => "I",
                "AAA" => "K",
                "AAG" => "K",
                "TTG" => "L",
                "TTA" => "L",
                "CTT" => "L",
                "CTC" => "L",
                "CTA" => "L",
                "CTG" => "L",
                "ATG" => "M",
                "AAT" => "N",
                "AAC" => "N",
                "CCT" => "P",
                "CCC" => "P",
                "CCA" => "P",
                "CCG" => "P",
                "CAA" => "Q",
                "CAG" => "Q",
                "CGT" => "R",
                "CGC" => "R",
                "CGA" => "R",
                "CGG" => "R",
                "AGA" => "R",
                "AGG" => "R",
                "TCT" => "S",
                "TCC" => "S",
                "TCA" => "S",
                "TCG" => "S",
                "AGT" => "S",
                "AGC" => "S",
                "ACT" => "T",
                "ACC" => "T",
                "ACA" => "T",
                "ACG" => "T",
                "GTT" => "V",
                "GTC" => "V",
                "GTA" => "V",
                "GTG" => "V",
                "TGG" => "W",
                "TAT" => "Y",
                "TAC" => "Y",
                "TAA" => "*",
                "TAG" => "*",
                "TGA" => "*"}
                self.geneticCode = gc
              end
              return self
            end

            def setPUDict	
              pu = { "GCT" => 0,
                "GCC" => 1,
                "GCA" => 0,
                "GCG" => 0,
                "TGT" => 0,
                "TGC" => 1,
                "GAT" => 0,
                "GAC" => 1,
                "GAA" => 0,
                "GAG" => 1,
                "TTT" => 0,
                "TTC" => 1,
                "GGT" => 0,
                "GGC" => 1,
                "GGA" => 0,
                "GGG" => 0,
                "CAT" => 0,
                "CAC" => 1,
                "ATT" => 0,
                "ATC" => 1,
                "ATA" => 0,
                "AAA" => 0,
                "AAG" => 1,
                "TTG" => 0,
                "TTA" => 0,
                "CTT" => 0,
                "CTC" => 1,
                "CTA" => 0,
                "CTG" => 1,
                "ATG" => 0,
                "AAT" => 0,
                "AAC" => 1,
                "CCT" => 0,
                "CCC" => 1,
                "CCA" => 0,
                "CCG" => 0,
                "CAA" => 0,
                "CAG" => 1,
                "CGT" => 1,
                "CGC" => 1,
                "CGA" => 0,
                "CGG" => 0,
                "AGA" => 0,
                "AGG" => 0,
                "TCT" => 0,
                "TCC" => 1,
                "TCA" => 0,
                "TCG" => 1,
                "AGT" => 0,
                "AGC" => 0,
                "ACT" => 0,
                "ACC" => 1,
                "ACA" => 0,
                "ACG" => 0,
                "GTT" => 0,
                "GTC" => 1,
                "GTA" => 0,
                "GTG" => 1,
                "TGG" => 0,
                "TAT" => 0,
                "TAC" => 1,
                "TAA" => 0,
                "TAG" => 0,
                "TGA" => 0 }
                self.puDict = pu
                return self
              end

              def degenCodonSite(aCodon, aSite)			#sets the synonymous "degeneracy" of each position in each codon; Li 1993
                h ={ 
                  "TAC" => [0,0,2],
                  "GTC" => [0,0,4],
                  "CTG" => [2,0,4],
                  "CAT" => [0,0,2],
                  "GCG" => [0,0,4],
                  "ACC" => [0,0,4],
                  "AGG" => [2,0,2],
                  "CCA" => [0,0,4],
                  "TTA" => [2,0,2],
                  "AAA" => [0,0,2],
                  "ATT" => [0,0,3],
                  "GGA" => [0,0,4],
                  "TGT" => [0,0,2],
                  "TCG" => [0,0,4],
                  "GAG" => [0,0,2],
                  "GCT" => [0,0,4],
                  "TGA" => [0,2,0],
                  "AGT" => [0,0,2],
                  "CGG" => [2,0,4],
                  "CAA" => [0,0,2],
                  "CCC" => [0,0,4],
                  "AAC" => [0,0,2],
                  "CTT" => [0,0,4],
                  "ATA" => [0,0,3],
                  "GGC" => [0,0,4],
                  "TTC" => [0,0,2],
                  "TAG" => [0,0,2],
                  "GTG" => [0,0,4],
                  "GAT" => [0,0,2],
                  "ACG" => [0,0,4],
                  "TCT" => [0,0,4],
                  "AGA" => [2,0,2],
                  "CGT" => [0,0,4],
                  "CTA" => [2,0,4],
                  "ATC" => [0,0,3],
                  "CAC" => [0,0,2],
                  "TGC" => [0,0,2],
                  "GCA" => [0,0,4],
                  "TAT" => [0,0,2],
                  "GTT" => [0,0,4],
                  "ACT" => [0,0,4],
                  "AGC" => [0,0,2],
                  "TCA" => [0,0,4],
                  "CGA" => [2,0,4],
                  "CCG" => [0,0,4],
                  "CTC" => [0,0,4],
                  "TTG" => [2,0,2],
                  "AAG" => [0,0,2],
                  "GGG" => [0,0,4],
                  "GAA" => [0,0,2],
                  "GCC" => [0,0,4],
                  "TAA" => [0,2,2],
                  "TGG" => [0,0,0],
                  "GTA" => [0,0,4],
                  "ACA" => [0,0,4],
                  "TCC" => [0,0,4],
                  "CGC" => [0,0,4],
                  "CAG" => [0,0,2],
                  "CCT" => [0,0,4],
                  "AAT" => [0,0,2],
                  "ATG" => [0,0,0],
                  "GGT" => [0,0,4],
                  "TTT" => [0,0,2],
                  "GAC" => [0,0,2]}

                  return h[aCodon][aSite]
                end


                def setSilentSiteDict      #again Nei and Gojorbori style creation, altered a little to accomodate stops; returns an array [sSites,nSites]
                  dict = Hash.new
                  codons = self.geneticCode.keys.sort
                  codons.delete("TAA")
                  codons.delete("TGA")
                  codons.delete("TAG") 
                  dict["TAA"] = [nil,nil]
                  dict["TAG"] = [nil,nil]
                  dict["TGA"] = [nil,nil]

                  #go through codons
                  codons.each{ | aCodon | 
                    s = 0
                    n = 9
                    origAA = self.geneticCode[aCodon]
                    0.upto(2){ | i |
                      temp = String.new
                      0.upto(2){  | j | temp[j,1] = aCodon[j,1]}
                      #set up possible states, delete actual
                      states = ["A","C","T","G"]
                      states.delete(temp[i,1])
                      states.each{ | aState |
                        temp[i,1] = aState
                        tempAA = self.geneticCode[temp]
                        #if the test aa is the same add a silent site
                        if tempAA == origAA
                          s += 1
                          n -= 1
                        else
                          if tempAA == '*'
                            n -= 1
                          end
                        end
                      }

                    }
                    dict[aCodon] = [s.to_f / 3, n.to_f / 3]  
                  }
                  self.silentSiteDict =  dict
                end


                def setCodonDists						#reads in and stores paths as Hash, this creates a dependancy on the file codonMatrix.txt
                  keys = self.geneticCode.keys.sort
                  keys.delete("TAA")
                  keys.delete("TGA")
                  keys.delete("TAG")
                  codonDictionary = Hash.new
                  keys.each{ | key | codonDictionary[key] = Hash.new }
                  aFile = File.new("/Users/adk/rubyStuff/codonMatrix.txt")  #/Network/Servers/i-dpgp.ucdavis.edu/Users/adk/rubyStuff/codonMatrix.txt")
                  aFile.each_line{ | line |
                    array = line.chomp.split("\t")
                    if  array.first != ""
                      from = array.first
                      array.each_with_index{ | item, index |
                        if index > 0
                          if item != "nil"
                            codonDictionary[from][keys[index - 1]] = item
                            codonDictionary[keys[index - 1]][from] = item
                          end
                        end
                      }
                    end
                  }
                  self.codonDists = codonDictionary
                end

                def degeneracyMatrix			#returns an array of arrays corresponding to the 'degeneracy' of each sequence. [nondegen., 2fold, 4fold, codonCount]
                  a = Array.new
                  self.codons.each{ | eachAlleleArray |
                    nd = 0
                    f2 = 0
                    f4 = 0
                    cc = 0
                    eachAlleleArray.each{ | codon |
                      if ! (codon.include?("N") or codon.include?("-"))
                        cc += 1
                        0.upto(2){ | i |
                          x = self.degenCodonSite(codon,i)
                          case x
                          when 0
                            nd += 1
                          when 2 
                            f2 += 1
                          when 3
                            f2 += 1
                          else
                            f4 += 1
                          end
                        }
                      end
                    }
                    a << [nd,f2,f4,cc]
                  }
                  return a
                end

                #silent site stuff, like sil poly
                #

                def gcThree    #average gc 3 values for seqs
                  freqs = Array.new
                  self.codons.each{ | eachRow |
                    gcCount = 0
                    codonCount = 0
                    eachRow.each{ |eachCodon |
                      if ! (eachCodon.include?("N") or eachCodon.include?("-"))
                        codonCount += 1
                        if eachCodon[2,1] =~ /[GC]/
                          gcCount += 1
                        end
                      end
                    }

                    if codonCount != 0
                      freqs << gcCount.to_f / codonCount
                    end
                  }


                  return freqs.sampleMean
                end

                def silentSites			#returns array  of silent sites
                  silSites = Array.new
                  self.codons.each{ | eachRow |
                    temp = Array.new
                    eachRow.each{ | eachCodon |
                      if ! ( eachCodon.include?("N") or eachCodon.include?("-") )
                        temp << self.silentSiteDict[eachCodon].first
                      end
                    }
                    silSites << temp.reject{ | x | x == nil }.sum
                  }
                  return silSites
                end

                def replacementSites
                  replSites = Array.new
                  self.codons.each{ | eachRow |
                    temp = Array.new
                    eachRow.each{ | eachCodon |
                      if ! ( eachCodon.include?("N") or eachCodon.include?("-") )
                        temp << self.silentSiteDict[eachCodon].last
                      end
                    }
                    replSites << temp.reject{ | x | x == nil}.sum
                  }
                  return replSites
                end

                def averageNumberSilentSites		#returns a float
                  return self.silentSites.sampleMean
                end

                def averageNumberReplacementSites
                  return self.replacementSites.sampleMean
                end

                def silentPi		#cleans ambiguites, uses parsimony counts and uses n&g silent site counts
                  sites = self.silReplMutations["silents"]
                  l = self.averageNumberSilentSites
                  oc = Array.new
                  sites.each{ | c |
                    siteArray = self.siteArrayClean(c)
                    siteSet = siteArray.uniq
                    ni = siteArray.size.to_f
                    if siteSet.size > 1
                      pSum = 0.0
                      siteSet.each{ | state |
                        p = (siteArray.occurrencesOf(state)).to_f/ ni
                        pSum += (p * p)
                      }
                      oc.push((1.0 - pSum) * (ni/(ni - 1))) 
                    end
                  }
                  if oc.empty?
                    return "0.0"
                  else
                    return (oc.sum/l)
                  end
                end

                def silReplPi2   #this routine uses the NG unweighted estimates of silent changes
                  comps = 0
                  silDifs = 0.0
                  replDifs = 0.0
                  #do all pairwise comparisons
                  0.upto(self.sampleSize - 2){ | i |
                    (i+1).upto(self.sampleSize - 1){ | j |
                      difs = self.silReplDifsAllPaths2Seqs(i, j)
                      if ! difs.include?(nil)
                        replDifs += difs[0]
                        silDifs += difs[1]
                        comps += 1
                      end
                    }
                  }
                  silDifs = silDifs / comps
                  if silDifs.nan?
                    silDifs = 0
                  end
                  replDifs = replDifs / comps
                  if replDifs.nan?
                    replDifs = 0
                  end
                  return [replDifs, silDifs]
                end


                def replacementPi					#cleans ambiguites and uses n&g silent site counts
                  sites = self.silReplMutations["replacements"]
                  l = self.averageNumberReplacementSites
                  oc = Array.new
                  sites.each{ | c |
                    siteArray = self.siteArrayClean(c)
                    siteSet = siteArray.uniq
                    ni = siteArray.size.to_f
                    if siteSet.size > 1
                      pSum = 0.0
                      siteSet.each{ | state |
                        p = (siteArray.occurrencesOf(state)).to_f/ ni
                        pSum += (p * p)
                      }
                      oc.push((1.0 - pSum) * (ni/(ni - 1))) 
                    end
                  }
                  if oc.empty?
                    return "0.0"
                  else
                    return (oc.sum/l)
                  end
                end

                def replacementTheta
                  sites = self.silReplMutations["replacements"]
                  l = self.averageNumberReplacementSites
                  s = sites.size / l
                  n = self.sampleSize
                  if n < 2
                    return 0.0
                  else
                    sum = 0.0
                    count = 1
                    while count < n
                      sum += (1.0/count)
                      count += 1
                    end
                  end
                  return s/sum
                end

                def silentTheta
                  sites = self.silReplMutations["silents"]
                  l = self.averageNumberSilentSites
                  s = sites.size / l
                  n = self.sampleSize
                  if n < 2
                    return 0.0
                  else
                    sum = 0.0
                    count = 1
                    while count < n
                      sum += (1.0/count)
                      count += 1
                    end
                  end
                  return s/sum
                end	

                def silentTajD

                  s = self.silReplMutations["silents"].size
                  if s < 3
                    return -666
                  else
                    n = self.sampleSize.to_f
                    if n < 4
                      return -666
                    else
                      k = self.silReplPi2[1] * self.averageNumberSilentSites
                      a1 = 0.0
                      a2 = 0.0
                      1.upto(n - 1){ | i |
                        a1 +=  (1.0/i)
                        a2 += (1.0/ (i*i))
                      }
                      b1 = (n + 1.0)/(3.0 * ( n - 1.0))
                      b2 = (2*( (n * n) + n + 3.0))/ (9.0 * n * (n - 1.0))
                      c1 = b1 - (1.0 / a1)
                      c2 = b2 - ((n + 2) / (a1 * n)) + (a2 / (a1 * a1))
                      e1 = c1 / a1
                      e2 = c2 / ((a1 * a1) + a2)
                      d = k - (s.to_f / a1)
                      stdev = Math.sqrt((e1 * s) + ((e2 * s) * (s - 1)))
                      return d/stdev
                    end
                  end
                end
                def replTajD

                  s = self.silReplMutations["replacements"].size
                  if s < 3
                    return -666
                  else
                    n = self.sampleSize.to_f
                    if n < 4
                      return -666
                    else
                      k = self.silReplPi2[0] * self.averageNumberReplacementSites
                      a1 = 0.0
                      a2 = 0.0
                      1.upto(n - 1){ | i |
                        a1 +=  (1.0/i)
                        a2 += (1.0/ (i*i))
                      }
                      b1 = (n + 1.0)/(3.0 * ( n - 1.0))
                      b2 = (2*( (n * n) + n + 3.0))/ (9.0 * n * (n - 1.0))
                      c1 = b1 - (1.0 / a1)
                      c2 = b2 - ((n + 2) / (a1 * n)) + (a2 / (a1 * a1))
                      e1 = c1 / a1
                      e2 = c2 / ((a1 * a1) + a2)
                      d = k - (s.to_f / a1)
                      stdev = Math.sqrt((e1 * s) + ((e2 * s) * (s - 1)))
                      return d/stdev
                    end
                  end
                end






              end






              class Array
                include ArrayADK
                include Enumerable
              end
              class SequenceMatrix
                include CodingSequence
              end

              class MultipleMatrix
                attr_reader :samples
                attr_writer :samples

                def initialize(aFile)
                  fileLines = File.readlines(aFile)
                  self.samples = Array.new
                  while ! fileLines.empty?
                    tempArray = fileLines.slice!(0,fileLines.index("//\n") + 1)
                    tempArray.pop
                    self.samples << SequenceMatrix.new.initializeFromString(tempArray)
                  end
                end
              end

              class SequenceMatrix
                attr_reader :matrix, :nameVector, :codingRegions, :readingFrame, :geneticCode, :codons, :features,:filenameID, :codonDists, :puDict, :silentSiteDict
                attr_writer :matrix, :nameVector, :codingRegions, :readingFrame, :geneticCode, :codons, :features, :filenameID, :codonDists, :puDict, :silentSiteDict


                def initializeFromFasta(aFilename)
                  aFile = File.new(aFilename)
                  array = Array.new
                  string = String.new("")
                  nameVector = Array.new
                  aFile.each_line{ | line |
                    if ! line.empty?
                      if line =~ /^>.+/
                        line.chomp!.slice!(/>/)
                        nameVector.push(line)
                        if string != ""
                          array.push(string)
                        end
                        string = String.new("")
                      else 
                        string << line.chomp.gsub(/ /,"").upcase
                      end
                    end }
                    array.push(string.upcase)
                    self.matrix = array
                    self.nameVector = nameVector
                    self.filenameID= aFilename
                    aFile.close
                    return self
                  end

                  def initializeFromString(aString)
                    if aString.class == String
                      lines = aString.split("\n")
                    else
                      lines = aString
                    end
                    array = Array.new
                    string = String.new("")
                    nameVector = Array.new
                    lines.each{ | line |
                      if ! line.empty?
                        if line =~ /^>.+/
                          line.slice!(/>/)
                          nameVector.push(line.chomp)
                          if string != ""
                            array.push(string)
                          end
                          string = String.new("")
                        else 
                          string << line.chomp.gsub(/ /,"").upcase
                        end
                      end }
                      array.push(string.upcase)
                      self.matrix = array
                      self.nameVector = nameVector
                      return self
                    end

                    def asCodingSequence(anArray,frame)
                      if frame == nil
                        self.readingFrame = 1
                      else
                        self.readingFrame = frame
                      end
                      if anArray != nil
                        self.codingRegions = anArray
                      else
                        self.codingRegions = [1,self.matrix[0].length]
                      end
                      self.setGeneticCode("standard")
                      self.setCodonDists
                      self.setPUDict
                      self.setSilentSiteDict
                      self.setCodons
                    end

                    #	def initializeFromGenbank(aFilename)
                    #		
                    #		aFile = File.new(aFilename)
                    #		array = Array.new
                    #		string = String.new("")
                    #		nameVector = Array.new
                    #		self.codingRegions = Array.new
                    #		featureDict = Hash.new
                    #		aFile.each_line{ | line |
                    #				if ! line.empty?
                    #					if line.strip =~ /^CDS/
                    #						line.chomp!.gsub(/d+/){ | match | self.codingRegions.push(match)
                    #															print self.codingRegions}
                    #					end	
                    #					if line.strip =~ /^ORIGIN/
                    #							array.push(string.upcase)
                    #						end
                    #						string = String.new("")
                    #					else 
                    #						string << line.chomp.gsub(" ","")
                    #					end
                    #				 
                    #				}
                    #		array.push(string.upcase)
                    #		self.matrix = array
                    #		self.nameVector = nameVector

                    #		self.filenameID= aFilename
                    #		self.readingFrame= 1
                    #		aFile.close
                    #		return self
                    #	end


                    def initializeWith(nameVector, matrix)
                      self.matrix = matrix
                      self.nameVector = nameVector
                      self.codingRegions = Array.new
                      return self
                    end

                    # manipulations
                    def concatenate			#returns a new instance of a SequenceMatrix with all of the sequences of the original concatenated into one long sequence
                      temp = String.new("")
                      name = self.nameVector[0]
                      self.matrix.each{ | eachString |
                        temp << eachString
                      }
                      return SequenceMatrix.new.initializeWith([name],[temp])
                    end

                    def sequenceSubSet(anArray)			#takes an array of sites, returns a new instance of a SequenceMatrix with only those sites
                      newMat = Array.new
                      0.upto(self.sampleSize - 1){ | r |
                        newSeq = String.new("")
                        anArray.each{ | c |
                          if c < self.length - 1
                            newSeq << self.matrix[r][c,1]
                          else
                            newSeq << "N"
                          end
                        }
                        newMat << newSeq
                      }
                      return SequenceMatrix.new.initializeWith(self.nameVector,newMat)
                    end

                    def sequenceMatrixSplit(anArray) #returns a new instance of a SequenceMatrix with only those seqs specified in array
                      newMat = Array.new
                      newNames = Array.new
                      anArray.each{ | anIndex |
                        if self.nameVector[anIndex] != nil
                          newNames << self.nameVector[anIndex].clone
                          newMat << self.matrix[anIndex].clone
                        end
                      }
                      return(SequenceMatrix.new.initializeWith(newNames,newMat.flatten))  
                    end

                    def	sequenceSlice(anArray)			#meant to take a start and stop array
                      newMat = Array.new
                      0.upto(self.sampleSize - 1){ | r |
                        l = anArray.last - anArray.first 
                        newSeq = self.matrix[r][anArray.first, l]
                        newMat << newSeq
                      }
                      return SequenceMatrix.new.initializeWith(self.nameVector,newMat)
                    end

                    def returnSequences(anArray)  #returns a new instance of sequenceMatrix with the sequences indicated in the array
                      newMat = Array.new
                      newNames = Array.new
                      anArray.each{ | x | 
                        newMat << self.matrix[x]
                        newNames << self.nameVector[x]
                      }
                      return SequenceMatrix.new.initializeWith(newNames,newMat)
                    end

                    def sampleSize
                      return @matrix.size
                    end

                    def length
                      return self.matrix[0].size
                    end

                    def siteSet(anIndex)
                      set = Array.new
                      @matrix.each{ | aSeq |
                        set.push(aSeq[anIndex,1])
                      }
                      set.uniq!
                      return set
                    end
                    def siteSetClean(anIndex)		#cleans out ambiguous bases
                      set = Array.new
                      @matrix.each{ | aSeq |
                        set.push(aSeq[anIndex,1])
                      }
                      set.uniq!
                      set.delete("N")
                      set.delete("-")
                      return set
                    end

                    def siteArray(anIndex)
                      array = Array.new
                      @matrix.each{ | aSeq |
                        array.push(aSeq[anIndex,1])
                      }
                      return array
                    end

                    def siteArrayClean(anIndex)
                      array = Array.new
                      @matrix.each{ | aSeq |
                        array.push(aSeq[anIndex,1])
                      }
                      array.delete("N")
                      array.delete("-")
                      return array
                    end

                    def aminoSet(anIndex)
                      set = Array.new
                      self.codons.each{ | anAllele |
                        if anAllele[anIndex] =~ /-/ 
                          set.push("-")
                        else if anAllele[anIndex] =~ /N/
                          set.push("N")
                        else set.push(self.geneticCode[anAllele[anIndex]])
                        end
                      end
                    }
                    return set.uniq
                  end

                  def aminoSetClean(anIndex)
                    set = Array.new
                    self.codons.each{ | anAllele |
                      if ! (anAllele[anIndex] =~ /-/ or anAllele[anIndex] =~ /N/)
                        set.push(self.geneticCode[anAllele[anIndex]])
                      end
                    }
                    return set.uniq
                  end

                  def translate  # returns a new instance of sequence matrix
                    self.setCodons
                    trans = SequenceMatrix.new
                    oc = Array.new
                    self.codons.each { | anAllele |	
                      temp = String.new("")
                      anAllele.each{ | aCodon |
                        if aCodon.include?("-") or aCodon.include?("N")
                          temp << "-"
                        else 
                          if self.geneticCode[aCodon.upcase]
                            temp << (self.geneticCode[aCodon.upcase])
                          else 
                            print aCodon
                          end
                        end }
                        oc.push(temp)
                      }
                      oc.each{ | row |
                        last = row.length - 1
                        if row[last,1] == "-" or row[last,1]  == "*"
                          row.slice!(last)
                        end
                      }
                      trans.matrix = oc
                      trans.nameVector = self.nameVector
                      return trans
                    end

                    def has_stops?    #data check. translates and looks for stop codons
                      trans = self.translate
                      count = 0
                      trans.matrix.each{ | r |
                        count += r.count("*") }
                        return count
                      end

                      def percentGaps			#this is for data checking, returns an array where each element is for an allele

                        array = Array.new
                        matrix.each{ | r |
                          gaps = r.count("-")
                          count = r.size.to_f
                          array.push(gaps/count)
                        }
                        return array
                      end

                      def percentAlignment   #returns a float with the percentage of nucleotides aligned in all alleles
                        count = 0
                        length = matrix[0].length
                        0.upto(length){ | i |
                          if self.siteArray(i).include?("-")
                            count += 1
                          end
                        }
                        return 1.0 - (count.to_f / length.to_f)
                      end

                      def averageSampleSize   #returns a float with the average sample size - N's thrown away
                        count = 0.0
                        n = self.sampleSize
                        length = self.length
                        0.upto(length - 1){ | i |
                          count += self.siteArrayClean(i).size
                        }
                        return count.to_f / length
                      end	

                      def distancesBetweenSegSites  	#returns an array of distances, caluculated right to left
                        dists = Array.new
                        oc = self.segSitesLocations
                        oc.reverse!
                        top = oc.size
                        count = 0
                        while count < top - 1
                          dists.push((oc[count].to_f - oc[count + 1].to_f))
                          count += 1
                        end
                        return dists
                      end

                      def getHaplotypes(aFlag) 		#returns an instance of sequenceMatrix with only the segregating sites included
                        if aFlag.nil?
                          locs = self.segSitesLocations
                        else 
                          locs = self.segSitesLocationsWithGaps
                        end
                        mat = Array.new
                        0.upto(self.sampleSize - 1){ | i | mat.push(String.new("")) }
                        locs.each{ | c |
                          0.upto(self.sampleSize - 1){ | r |
                            mat[r]<<(self.matrix[r])[c,1]
                          }
                        }
                        newSeqMat = SequenceMatrix.new
                        newSeqMat.nameVector = self.nameVector
                        newSeqMat.matrix = mat
                        newSeqMat.readingFrame = 1
                        return newSeqMat
                      end

                      def getHaplotypesPretty(aFlag)	# returns an instance of sequenceMatrix with only the segregating sites but with .'s for the same as ref

                        haps = self.getHaplotypes(aFlag)
                        ref = haps.matrix[0]
                        mat = Array.new			
                        0.upto(haps.sampleSize - 1){ | i | mat.push(String.new("")) }
                        mat[0] = ref
                        0.upto(ref.length - 1){ | c |
                          1.upto(haps.sampleSize - 1){ | r |
                            if (haps.matrix[r])[c,1] == ref[c,1]
                              mat[r]<< "."
                            else
                              mat[r]<<(haps.matrix[r])[c,1]
                            end
                          }
                        }
                        newSeqMat = SequenceMatrix.new
                        newSeqMat.nameVector = haps.nameVector
                        newSeqMat.matrix = mat
                        newSeqMat.readingFrame = 1
                        return newSeqMat
                      end

                      def majorAlleleFreqSite(aSite)
                        tempArray = self.siteArrayClean(aSite)
                        set = self.siteSetClean(aSite)
                        n = tempArray.size.to_f
                        major = 0.0
                        set.each{ | each | 
                          if tempArray.occurrencesOf(each) > major 
                            major = tempArray.occurrencesOf(each)
                          end
                        }
                        return major/n
                      end


                      def majorAlleleFreqSiteDirty(aSite)
                        tempArray = self.siteArray(aSite)
                        set = self.siteSet(aSite)
                        n = tempArray.size.to_f
                        major = 0.0
                        set.each{ | each | 
                          if tempArray.occurrencesOf(each) > major 
                            major = tempArray.occurrencesOf(each)
                          end
                        }
                        return major/n
                      end

                      def minorAlleleFreqSite(aSite)
                        tempArray = self.siteArrayClean(aSite)
                        set = self.siteSetClean(aSite)
                        n = tempArray.size.to_f
                        minor = self.sampleSize
                        set.each{ | each | 
                          if tempArray.occurrencesOf(each) < minor 
                            minor = tempArray.occurrencesOf(each)
                          end
                        }
                        return minor
                      end
                      def minorAlleleFreqSiteDirty(aSite)
                        tempArray = self.siteArray(aSite)
                        set = self.siteSet(aSite)
                        n = tempArray.size.to_f
                        minor = self.sampleSize
                        set.each{ | each | 
                          if tempArray.occurrencesOf(each) < minor 
                            minor = tempArray.occurrencesOf(each)
                          end
                        }
                        return minor/n
                      end

                      def majorAlleleNumberSite(aSite)
                        tempArray = self.siteArray(aSite)
                        set = self.siteSet(aSite)
                        edit = tempArray.delete_if { |x | x =="N" }
                        major = 0.0
                        set.each{ | each | 
                          if edit.occurrencesOf(each) > major 
                            major = edit.occurrencesOf(each)
                          end
                        }
                        return major
                      end


                      def majorAlleleStateSite(aSite)
                        tempArray = self.siteArray(aSite)
                        set = self.siteSet(aSite)
                        edit = tempArray.delete_if { |x | x =="N" }
                        n = edit.size.to_f
                        major = 0
                        state = nil
                        set.each{ | each | 
                          if edit.occurrencesOf(each) > major 
                            major = edit.occurrencesOf(each)
                            state = each
                          end
                        }
                        return state
                      end

                      def stateNumberSite(aSite, aState)
                        tempArray = self.siteArray(aSite)
                        return tempArray.occurrencesOf(aState)
                      end

                      def rSquaredSites(aSite1, aSite2)

                        pA = self.majorAlleleFreqSite(aSite1)
                        stateA = self.majorAlleleStateSite(aSite1)
                        pB = self.majorAlleleFreqSite(aSite2)
                        stateB = self.majorAlleleStateSite(aSite2)

                        array1 = self.siteArrayClean(aSite1)
                        array2 = self.siteArrayClean(aSite2)
                        count = 0.0
                        siteCount = [array1.length,array2.length].min

                        0.upto(self.sampleSize - 1) { | i |
                          if (array1[i] == stateA) and (array2[i] == stateB)
                            count += 1.0
                          end
                        }
                        pAB = count/siteCount
                        d = pAB - (pA * pB)
                        denom = Math.sqrt(pA * (1 - pA) * pB * (1 - pB))
                        r = d / denom
                        return (r * r)
                      end

                      def reverseComplement   #returns sequence matrix

                        anArray = Array.new
                        self.matrix.each{ | string |
                          temp = string.reverse.tr('ATGC','TACG')
                          anArray.push(temp)
                        }
                        return SequenceMatrix.new.initializeWith(self.nameVector, anArray)
                      end

                      def relativePosition(aSequenceNumber, anIndex)			#returns the position of anIndex relative to the gaps inserted during alignment in sequence aSequenceNumber, useful for annotations but note it is zero based
                        tempIndex = relIndex = 0
                        string = self.matrix[aSequenceNumber]
                        while tempIndex <= anIndex
                          if string[relIndex,1] == "-"
                            relIndex += 1
                          else
                            relIndex += 1
                            tempIndex += 1
                          end

                        end
                        return relIndex - 1
                      end

                      def padAlignment		#evens sequence length with N's from the back
                        max = 0
                        self.matrix.each{ | eachString |
                          if eachString.length > max
                            max = eachString.length
                          end
                        }
                        self.matrix.each{ | eachString |
                          if eachString.length < max
                            l = max - eachString.length
                            l.times{ | i | eachString << "N" }
                          end
                        }
                      end

                      def padAlignmentLength(anInt)		#adds sequence length with N's from the back to length anInt
                        self.matrix.each{ | eachString |
                          if eachString.length < anInt
                            l = anInt - eachString.length
                            l.times{ | i | eachString << "N" }
                          end
                        }
                      end

                      # conversions
                      def codonsToPaml		#this attempts to cut off termination codons- still needs debugging
                        fn = self.filenameID, ".paml"
                        ws = File.new(fn.to_s, "w")
                        n = self.sampleSize 
                        size = self.codons[0].size - 1
                        l = (size * 3) + 3
                        if (self.aminoSet((size - 1))).include?("*")
                          l = l - 3
                          ws.print(n.to_s,"\t",l.to_s,"\n")
                          self.codons.each_index{ | i |
                            name = nameVector[i].to_s
                            ws.print(name.slice(0,10),"  ")
                            0.upto(self.codons[0].size - 2){ | c |
                              ws.print(self.codons[i][c])
                            }
                            ws.print("\n") }
                          else
                            ws.print(n.to_s,"\t",l.to_s,"\n")
                            self.codons.each_index{ | i |
                              name = nameVector[i].to_s
                              ws.print(name.slice(0,10),"  ")
                              0.upto(self.codons[0].size - 1){ | c |
                                ws.print(self.codons[i][c])
                              }
                              ws.print("\n") }
                            end
                            ws.close
                          end

                          def sequenceNameHash		#returns a hash with the sequence names as keys and the sequence strings as values
                            i = 0
                            aHash = Hash.new
                            self.nameVector.each{ | key |
                              aHash[key] = self.matrix[i]
                              i += 1
                            }
                            return aHash
                          end

                          def returnFasta    #spits out a new fasta file reflecting any changes that might have occured
                            fn = self.filenameID, ".new"
                            ws = File.new(fn.to_s, "w")
                            self.nameVector.each_index{ | i |
                              ws.print(">",nameVector[i],"\n",matrix[i],"\n")
                            }
                            ws.close
                          end

                          def returnFastaNamed(aName)    #spits out a new fasta file reflecting any changes that might have occured
                            ws = File.new(aName, "w")
                            self.nameVector.each_index{ | i |
                              ws.print(">",nameVector[i],"\n",matrix[i],"\n")
                            }
                            ws.close
                          end

                          def outputFasta		
                            i = 0
                            self.nameVector.each{ | name |
                              print ">",name,"\n"
                              count = 60
                              0.upto(self.matrix[i].length - 1){ | j |
                                print self.matrix[i][j,1]
                                count -= 1
                                if count == 0
                                  print "\n"
                                  count = 60
                                end
                              }
                              i += 1
                              print "\n\n"
                            }
                          end


                          # polymorphism stats

                          def segSites   #this ignores N's but keeps -'s. will be different from estimator S's
                            s = 0
                            count = 0
                            while count < self.matrix[0].size  
                              temp = self.siteSet(count)
                              temp.delete("N")
                              if 1 < temp.size
                                s += 1
                              end
                              count += 1
                            end
                            return s
                          end

                          def segSitesKillN
                            s = 0
                            count = 0
                            while count < self.matrix[0].size  
                              temp = self.siteSet(count)
                              if temp.include?("N") or temp.include?("-")
                                count += 1
                              else
                                if 1 < temp.size
                                  s += 1
                                end
                                count += 1
                              end
                            end
                            return s
                          end

                          def sitesNArray   #returns array length n of number of sites at each sampleSize
                            count = 0
                            n = self.sampleSize
                            oc = Array.new
                            n.times{ | i | oc << 0 }
                            while count < self.matrix[0].size  
                              ni = self.siteArrayClean(count).size
                              oc[ni - 1] += 1
                              count += 1
                            end
                            return oc
                          end    

                          def segSitesNArray   #returns array length N of number of segSites at each sampleSize
                            count = 0
                            n = self.sampleSize
                            oc = Array.new
                            n.times{ | i | oc << 0 }
                            while count < self.matrix[0].size  
                              temp = self.siteSetClean(count)
                              if 1 < temp.size
                                ni = self.siteArrayClean(count).size
                                oc[ni - 1] += 1
                              end
                              count += 1
                            end
                            return oc
                          end

                          def segSitesNArrayWindow(start, fin)   #returns array length N of number of segSites at each sampleSize, note zero indexing- doesn't look at seq at position "fin"
                            count = start
                            n = self.sampleSize
                            oc = Array.new
                            n.times{ | i | oc << 0 }
                            while count < fin  
                              temp = self.siteSetClean(count)
                              if 1 < temp.size
                                ni = self.siteArrayClean(count).size
                                oc[ni - 1] += 1
                              end
                              count += 1
                            end
                            return oc
                          end

                          def segSitesLocations 			#returns array of locations from segSitesKillN zero indexed
                            oc = Array.new				# note that this is for sites counted in estimators like theta 
                            count = 0
                            while count < self.matrix[0].size  
                              temp = self.siteSet(count)
                              if temp.include?("N") or temp.include?("-")
                                count += 1
                              else
                                if 1 < temp.size
                                  oc.push(count)
                                end
                                count += 1
                              end
                            end
                            return oc
                          end

                          def segSitesLocationsWithGaps 			#returns array of locations of segSites with gaps includes (still no N's)
                            oc = Array.new
                            count = 0
                            while count < self.matrix[0].size  
                              temp = self.siteSet(count)
                              if temp.include?("N")
                                count += 1
                              else
                                if 1 < temp.size
                                  oc.push(count)
                                end
                                count += 1
                              end
                            end
                            return oc
                          end

                          def segSitesLocationsAllSites 			#returns array of locations from segSitesKillN zero indexed
                            oc = Array.new				# note that this is for all columns (even those that contain -'s or N's)
                            count = 0
                            while count < self.matrix[0].size  
                              temp = self.siteSetClean(count)
                              if 1 < temp.size
                                oc.push(count)
                              end
                              count += 1			
                            end
                            return oc
                          end

                          def segSitesLocationsAllSitesNoSingletons 	#returns array of locations from segSitesKillN zero indexed
                            oc = Array.new				# note that this is for all columns (even those that contain -'s or N's)
                            count = 0
                            while count < self.matrix[0].size  
                              temp = self.siteSetClean(count)
                              if 1 < temp.size and minorAlleleFreqSite(count) > 1
                                oc.push(count)
                              end
                              count += 1			
                            end
                            return oc
                          end

                          def segSitesLocationsAllSitesWithGaps 			# same as above but count gaps. returns array of locations from segSitesKillN zero indexed
                            oc = Array.new								# note that this is for all columns (even those that contain -'s or N's)
                            count = 0
                            while count < self.matrix[0].size  
                              temp = self.siteSet(count)
                              temp.delete("N")
                              if 1 < temp.size
                                oc.push(count)
                              end
                              count += 1

                            end
                            return oc
                          end

                          def fractionSegSitesKillN		#used in calculated Watterson's theta
                            s = 0
                            count = 0
                            siteCount = 0
                            while count < self.matrix[0].size  
                              temp = self.siteSet(count)
                              if temp.include?("N") or temp.include?("-")
                                count += 1
                              else
                                siteCount += 1
                                if 1 < temp.size
                                  s += 1
                                end
                                count += 1
                              end
                            end
                            return s.to_f/siteCount
                          end

                          def SequenceMatrix.a1(n)			#this is the harmonic sum in Watterson's estimator
                            sum = 0
                            1.upto(n - 1){ | i | sum += 1.0/i }
                            return sum
                          end

                          def thetaKillN 		#this is Watterson's estimator where sites which contain indels and N's are excluded
                            s = self.fractionSegSitesKillN
                            n = self.sampleSize
                            if n < 2
                              return 0.0
                            else
                              sum = 0.0
                              count = 1
                              while count < n
                                sum += (1.0/count)
                                count += 1
                              end
                            end
                            return s/sum
                          end

                          def thetaMissingData		#adapted version of Watterson's estimator
                            sArray = self.segSitesNArray
                            l = self.length
                            sum = 0
                            1.upto(sArray.size - 1){ | i | sum += sArray[i].to_f / SequenceMatrix.a1(i + 1) }
                            return sum / l	
                          end	

                          def thetaMissingDataWindow(start, fin)
                            sArray = self.segSitesNArrayWindow(start, fin)
                            l = fin - start
                            sum = 0
                            1.upto(sArray.size - 1){ | i | sum += sArray[i].to_f / SequenceMatrix.a1(i + 1) }
                            return sum / l	
                          end	

                          def piKillN  #this is Tajima's (1983) estimator 
                            n = self.sampleSize
                            oc = Array.new
                            l = self.matrix[0].size
                            index = 0
                            while index < self.matrix[0].size
                              siteArray = self.siteArray(index)
                              if siteArray.include?("-") or siteArray.include?("N")
                                l -= 1
                                index += 1
                              else
                                siteSet = siteSet(index)
                                ni = siteArray.size.to_f
                                if siteSet.size > 1
                                  pSum = 0.0
                                  siteSet.each{ | state |
                                    p = (siteArray.occurrencesOf(state)).to_f/ ni
                                    pSum += (p * p)
                                  }
                                  oc.push((1.0 - pSum) * (ni/(ni - 1))) 
                                end
                                index += 1
                              end
                            end
                            return (oc.sum/l)
                          end

                          def piClean				 #pi where sample size varies per column, also cleans gaps
                            oc = Array.new
                            l = self.matrix[0].size
                            index = 0
                            while index < self.matrix[0].size
                              siteArray = self.siteArrayClean(index)
                              siteSet = siteArray.uniq
                              ni = siteArray.size.to_f
                              if siteSet.size > 1
                                pSum = 0.0
                                siteSet.each{ | state |
                                  p = (siteArray.occurrencesOf(state)).to_f/ ni
                                  pSum += (p * p)
                                }
                                oc.push((1.0 - pSum) * (ni/(ni - 1))) 
                              end
                              index += 1

                            end
                            if oc.empty?
                              return "0.0"
                            else
                              return (oc.sum/l)
                            end
                          end

                          def piCleanNoSingletons				 #pi where sample size varies per column, also cleans gaps
                            oc = Array.new
                            l = self.matrix[0].size
                            self.segSitesLocationsAllSitesNoSingletons.each{ | index |
                              siteArray = self.siteArrayClean(index)
                              siteSet = siteArray.uniq
                              ni = siteArray.size.to_f
                              if siteSet.size > 1
                                pSum = 0.0
                                siteSet.each{ | state |
                                  p = (siteArray.occurrencesOf(state)).to_f/ ni
                                  pSum += (p * p)
                                }
                                oc.push((1.0 - pSum) * (ni/(ni - 1))) 
                              end
                            }
                            if oc.empty?
                              return "0.0"
                            else
                              return (oc.sum/l)
                            end
                          end

                          def piCleanWindow(start, fin)				 #pi where sample size varies per column, also cleans gaps, zero indexing- doesn't look at site "fin"
                            oc = Array.new
                            l = fin - start 
                            index = start
                            nCount = 0
                            while index < fin
                              siteArray = self.siteArrayClean(index)
                              if siteArray.empty?
                                nCount += 1
                              else
                                siteSet = siteArray.uniq
                                ni = siteArray.size.to_f
                                if siteSet.size > 1
                                  pSum = 0.0
                                  siteSet.each{ | state |
                                    p = (siteArray.occurrencesOf(state)).to_f/ ni
                                    pSum += (p * p)
                                  }
                                  oc.push((1.0 - pSum) * (ni/(ni - 1))) 
                                end
                              end
                              index += 1
                            end
                            if nCount > l * 0.5
                              return "NA"
                            else
                              return (oc.sum/l)
                            end
                          end
                          def tajimasK 	#this is pi total for Tajima's D
                            n = self.sampleSize
                            oc = Array.new
                            l = self.matrix[0].size
                            index = 0
                            while index < self.matrix[0].size
                              siteArray = self.siteArray(index)
                              if siteArray.include?("-") or siteArray.include?("N")
                                index += 1
                              else
                                siteSet = siteSet(index)
                                ni = siteArray.size.to_f
                                if siteSet.size > 1
                                  pSum = 0.0
                                  siteSet.each{ | state |
                                    p = (siteArray.occurrencesOf(state)).to_f/ ni
                                    pSum += (p * p)
                                  }
                                  oc.push((1.0 - pSum) * (ni/(ni - 1))) 
                                end
                                index += 1
                              end
                            end
                            return oc.sum
                          end

                          def thetaH(state) 	#Fay and Wu's H estimator, variable is the derived state
                            n = self.sampleSize
                            l = self.matrix[0].size
                            index = 0
                            pSum = 0.0
                            while index < self.matrix[0].size
                              siteArray = self.siteArray(index)
                              if siteArray.include?("-") or siteArray.include?("N")
                                index += 1
                              else
                                siteSet = siteSet(index)
                                ni = siteArray.size.to_f
                                if siteSet.size > 1
                                  p = (siteArray.occurrencesOf(state)).to_f
                                  pSum += (p * p)
                                end
                                index += 1
                              end
                            end
                            return (pSum * 2.0) / (n * (n-1.0))
                          end

                          def hFay(state) 	#Fay and Wu's H 
                            n = self.sampleSize
                            nnm1 = n / (n-1.0)
                            l = self.matrix[0].size
                            index = 0
                            pSum = 0.0
                            while index < self.matrix[0].size
                              siteArray = self.siteArray(index)
                              if siteArray.include?("-") or siteArray.include?("N")
                                index += 1
                              else
                                siteSet = siteSet(index)
                                ni = siteArray.size.to_f
                                if siteSet.size > 1
                                  p = (siteArray.occurrencesOf(state)).to_f / n
                                  pSum += 2.0 * p * (2.0 * p - 1) * nnm1
                                end
                                index += 1
                              end
                            end
                            return (-pSum)
                          end

                          def numberHaplotypes
                            haps = self.getHaplotypes(nil)
                            return haps.matrix.uniq.size
                          end

                          def haplotypeDiversity
                            haps = self.getHaplotypes(nil)
                            set = haps.matrix.uniq
                            sum = 0.0
                            set.each{ | string | 
                              p = haps.matrix.select{ | x | x == string }.size.to_f / haps.sampleSize.to_f
                              sum += p * p
                            }
                            return 1.0 - sum
                          end

                          #foldedSFS- returns array of site frequency spectrum assuming complete sample size
                          def sfs
                            s = self.segSitesLocationsAllSites
                            sfs = Array.new
                            0.upto(self.sampleSize - 1){ | i | sfs[i] = 0}
                            s.each{ | i |
                              sfs[self.minorAlleleFreqSite(i)] += 1
                            }
                            return sfs
                          end

                          def tajimasDKillN
                            s = self.segSitesKillN
                            if s < 3
                              return "na"
                            else
                              n = self.sampleSize.to_f
                              if n < 4
                                return "na"
                              else
                                a1 = 0.0
                                a2 = 0.0
                                1.upto(n - 1){ | i |
                                  a1 +=  (1.0/i)
                                  a2 += (1.0/ (i*i))
                                }
                                b1 = (n + 1.0)/(3.0 * ( n - 1.0))
                                b2 = (2*( (n * n) + n + 3.0))/ (9.0 * n * (n - 1.0))
                                c1 = b1 - (1.0 / a1)
                                c2 = b2 - ((n + 2) / (a1 * n)) + (a2 / (a1 * a1))
                                e1 = c1 / a1
                                e2 = c2 / ((a1 * a1) + a2)
                                d = self.tajimasK - (s.to_f / a1)
                                stdev = Math.sqrt((e1 * s) + ((e2 * s) * (s - 1)))
                                return d/stdev
                              end
                            end
                          end

                          def tajimasDHack
                            s = self.segSitesLocationsAllSites.size
                            if s < 3
                              return "na"
                            else
                              n = self.sampleSize.to_f
                              if n < 4
                                return "na"
                              else
                                a1 = 0.0
                                a2 = 0.0
                                1.upto(n - 1){ | i |
                                  a1 +=  (1.0/i)
                                  a2 += (1.0/ (i*i))
                                }
                                b1 = (n + 1.0)/(3.0 * ( n - 1.0))
                                b2 = (2*( (n * n) + n + 3.0))/ (9.0 * n * (n - 1.0))
                                c1 = b1 - (1.0 / a1)
                                c2 = b2 - ((n + 2) / (a1 * n)) + (a2 / (a1 * a1))
                                e1 = c1 / a1
                                e2 = c2 / ((a1 * a1) + a2)
                                d = (self.piClean * self.length) - (s.to_f / a1)
                                stdev = Math.sqrt((e1 * s) + ((e2 * s) * (s - 1)))
                                return d/stdev
                              end
                            end
                          end

                          def tajimasDHackNoSingletons
                            s = self.segSitesLocationsAllSitesNoSingletons.size
                            if s < 3
                              return "na"
                            else
                              n = self.sampleSize.to_f
                              if n < 4
                                return "na"
                              else
                                a1 = 0.0
                                a2 = 0.0
                                1.upto(n - 1){ | i |
                                  a1 +=  (1.0/i)
                                  a2 += (1.0/ (i*i))
                                }
                                b1 = (n + 1.0)/(3.0 * ( n - 1.0))
                                b2 = (2*( (n * n) + n + 3.0))/ (9.0 * n * (n - 1.0))
                                c1 = b1 - (1.0 / a1)
                                c2 = b2 - ((n + 2) / (a1 * n)) + (a2 / (a1 * a1))
                                e1 = c1 / a1
                                e2 = c2 / ((a1 * a1) + a2)
                                d = (self.piCleanNoSingletons * self.length) - (s.to_f / a1)
                                stdev = Math.sqrt((e1 * s) + ((e2 * s) * (s - 1)))
                                return d/stdev
                              end
                            end
                          end

                          def slidingWindowPi(windowSize, offset)    #note this uses piClean, prints array with [site, pi]
                            i = 0
                            l = self.length - windowSize
                            oc = Array.new
                            while i <= l
                              oc << [i, self.piCleanWindow(i, i + windowSize)]
                              i += offset
                            end
                            oc.reject!{ | each | each[1] == "NA" }
                            oc.each{ | each | print each[0],"\t",each[1],"\n"}
                          end

                          def slidingWindowTheta(windowSize, offset)    #note this uses thetaMissingData, prints array with [site, theta]
                            i = 0
                            l = self.length - windowSize
                            oc = Array.new
                            while i <= l
                              oc << [i, self.thetaMissingDataWindow(i, i + windowSize)]
                              i += offset
                            end
                            oc.each{ | each | print each[0],"\t",each[1],"\n"}
                          end

                          def polySummary    #returns an array with the name of the file, n, length, S, Stotal, H, thetaKillN, piKillN, tajimasDKillN
                            array = Array.new
                            array.push(self.filenameID)
                            array.push(self.sampleSize)
                            array.push(self.matrix[0].length)
                            array.push(self.segSitesKillN)
                            array.push(self.segSites)
                            array.push(self.numberHaplotypes)
                            array.push(self.haplotypeDiversity)
                            array.push(self.thetaKillN)
                            array.push(self.piKillN)
                            array.push(self.tajimasDKillN)
                            return array
                          end

                          def printPolySummary(anInt)
                            if anInt.nil?
                              sum = self.polySummary
                              print "file ID \tN \tbp \tS \tStotal \tH \thaploDiveristy \tthetaKillN \tpiKillN \ttajimasDKillN\n"
                              sum.each{| x | print x,"\t" }
                              print "\n"
                            else if anInt == "-p"
                              sum = self.polySummary
                              string = "ms ",sum[1].to_s," 10000 -s ",sum[3].to_s," | sample_stats > ~/rubyStuff/simTemp"
                              #	print "now running simulations: ",string,"\n"
                              system(string.to_s)
                              aFile = File.new("/Users/adk/rubyStuff/simTemp","r")
                              lines = aFile.readlines
                              array = Array.new
                              lines.each{| x | tokens = x.chomp.split 
                                array.push(tokens[5].to_f)}
                                array.sort!
                                sum.push(array[499])
                                sum.push(array[9499])
                                print "file ID \tN \tbp \tS \tStotal \tH \thaploDiveristy  \tthetaKillN \tpiKillN \ttajimasDKillN\tlower 95% d\tupper 95% d\n"
                                sum.each{| x | print x,"\t" }
                                print "\n"
                                `rm ~/rubyStuff/simTemp`
                              end
                            end
                          end	

                          #
                          # divergence (still not sure if this should be moved to separate class)
                          #

                          def averagePairwiseDifferencesPerSite			#assumes all seqs are the same length, kills "-"s and "N"s

                            difPerSite = Array.new
                            n = self.sampleSize - 1 
                            width = self.matrix[0].size - 1
                            0.upto(n - 1){ | i |
                              (i+1).upto(n){ | j |
                                siteCount = 0
                                diffs = 0
                                0.upto(width){ | c |
                                  testArray = [self.matrix[i][c,1],self.matrix[j][c,1]]
                                  if ! testArray.include?("-") or ! testArray.include?("N")
                                    siteCount += 1
                                    if testArray.uniq.size > 1
                                      diffs += 1
                                    end
                                  end
                                }
                                difPerSite << diffs.to_f / siteCount
                              }
                            }
                            return difPerSite.sampleMean
                          end

                          def averagePairwiseDifferencesArray			#assumes all seqs are the same length, kills "-"s and "N"s

                            difPerSite = Array.new
                            n = self.sampleSize - 1
                            width = self.matrix[0].size - 1
                            0.upto(n - 1){ | i |
                              (i+1).upto(n){ | j |
                                siteCount = 0
                                diffs = 0
                                0.upto(width){ | c |
                                  testArray = [self.matrix[i][c,1],self.matrix[j][c,1]]
                                  if ! testArray.include?("-") and ! testArray.include?("N")
                                    siteCount += 1
                                    if testArray.uniq.size > 1
                                      diffs += 1
                                    end
                                  end
                                }
                                difPerSite << diffs.to_f 
                              }
                            }
                            return difPerSite
                          end

                          def averagePairwiseDifferencesPerSiteSites(anArray)			#assumes all seqs are the same length, doesn't kill N's but kills "-"s

                            difPerSite = Array.new
                            n = self.sampleSize - 1 
                            0.upto(n - 1){ | i |
                              (i+1).upto(n){ | j |
                                siteCount = 0
                                diffs = 0
                                anArray.each{ | c |
                                  testArray = [self.matrix[i][c,1],self.matrix[j][c,1]]
                                  if ! testArray.include?("-")
                                    siteCount += 1
                                    if testArray.uniq.size > 1
                                      diffs += 1
                                    end
                                  end
                                }
                                if siteCount == 0
                                  difPerSite << "NA"
                                else
                                  difPerSite << diffs.to_f/siteCount
                                end
                              }
                            }
                            temp = difPerSite.checkArray
                            if temp
                              return temp.sampleMean
                            else
                              return "NA"
                            end
                          end

                          def averagePairwiseDifferencesSites(anArray)			#assumes all seqs are the same length, doesn't kill N's but kills "-"s

                            difPerSite = Array.new
                            n = self.sampleSize - 1 
                            0.upto(n - 1){ | i |
                              (i+1).upto(n){ | j |
                                diffs = 0
                                anArray.each{ | c |
                                  testArray = [self.matrix[i][c,1],self.matrix[j][c,1]]
                                  if ! testArray.include?("-") or ! testArray.include?("N")
                                    if testArray.uniq.size > 1
                                      diffs += 1
                                    end
                                  end
                                }
                                difPerSite << diffs.to_f
                              }
                            }
                            return difPerSite.sampleMean
                          end


                          #
                          # stats 
                          #
                          #

                          def pairwiseRSquaredByDistance		#returns an Array

                            temp = self.segSitesLocationsAllSitesWithGaps

                            sites = temp.delete_if{ | x | self.majorAlleleFreqSite(x) > 0.95 }
                            n = sites.size - 1
                            counter = 0
                            anArray = Array.new
                            i = 0
                            while i < n
                              j = i + 1
                              while j <= n
                                print "comparison ",counter.to_s,"\n"
                                counter+=1
                                temp = Array.new
                                temp.push(self.rSquaredSites(sites[i],sites[j]))
                                temp.push(sites[j] - sites[i])
                                anArray.push(temp)
                                j += 1
                              end
                              i += 1
                            end
                            return anArray
                          end

                          def pairwiseRSquared		#returns lower diag. matrix
                            temp = self.segSitesLocationsAllSitesWithGaps
                            sites = temp.delete_if{ | x | self.majorAlleleFreqSite(x) > 0.95 }
                            n = sites.size - 1
                            counter = 0
                            anArray = Array.new
                            i = 0
                            while i < n
                              j = i + 1
                              temp = Array.new(sites.length,"NA")
                              while j <= n    
                                temp[j] = self.rSquaredSites(sites[i],sites[j])
                                #temp.push(sites[j] - sites[i])
                                j += 1
                              end
                              anArray.push(temp)
                              i += 1
                            end
                            return anArray
                          end


                          ###############################################	
                          # Protein Sequence Stuff	
                          def segSitesLocationsAA 			#returns array of locations from segSitesKillN zero indexed
                            oc = Array.new				# note that this is for sites counted in estimators like theta 
                            count = 0
                            while count < self.matrix[0].size  
                              temp = self.siteSet(count)
                              if temp.include?("X") or temp.include?("-")
                                count += 1
                              else
                                if 1 < temp.size
                                  oc.push(count)
                                end
                                count += 1
                              end
                            end
                            return oc
                          end

                          def segSitesLocationsWithGapsAA 			#returns array of locations of segSites with gaps includes (still no N's)
                            oc = Array.new
                            count = 0
                            while count < self.matrix[0].size  
                              temp = self.siteSet(count)
                              if temp.include?("X")
                                count += 1
                              else
                                if 1 < temp.size
                                  oc.push(count)
                                end
                                count += 1
                              end
                            end
                            return oc
                          end	

                          def getHaplotypesAA(aFlag) 		#returns an instance of sequenceMatrix with only the segregating sites included
                            if aFlag.nil?
                              locs = self.segSitesLocationsAA
                            else 
                              locs = self.segSitesLocationsWithGapsAA
                            end
                            mat = Array.new
                            0.upto(self.sampleSize - 1){ | i | mat.push(String.new("")) }
                            locs.each{ | c |
                              0.upto(self.sampleSize - 1){ | r |
                                mat[r]<<(self.matrix[r])[c,1]
                              }
                            }
                            newSeqMat = SequenceMatrix.new
                            newSeqMat.nameVector = self.nameVector
                            newSeqMat.matrix = mat
                            newSeqMat.readingFrame = 1
                            return newSeqMat
                          end	

                          ###########################################
                          ###   SFS stuff
                          ##
                          ##unfoldedSFS-- for ms output
                          def unfoldedSFSHudson
                            s = self.segSitesLocationsAllSites
                            sfs = Array.new
                            0.upto(self.sampleSize - 1){ | i | sfs[i] = 0}
                            s.each{ | i |
                              sfs[self.stateNumberSite(i,"1")] += 1
                            }
                            return sfs
                          end

                          ##Uses Achaz's "system" to calculate theta_w based on sfs
                          def sfs2ThetaW(sfsArray)
                            weightSum = 0.0
                            n = sfsArray.length
                            1.upto(n-1){|i|
                              weightSum += (1.0 / i)
                            }
                            thetaW = 0.0
                            1.upto(n-1){|i|
                              thetaW += sfsArray[i] * i * (1.0 / i)
                            }
                            return(thetaW * (1.0 / weightSum))
                          end

                          ##Uses Achaz's "system" to calculate theta_w based on sfs
                          # doesn't include singletons
                          def sfs2ThetaWNoSingletons(sfsArray)
                            weightSum = 0.0
                            n = sfsArray.length
                            2.upto(n-1){|i|
                              weightSum += (1.0 / i)
                            }
                            thetaW = 0.0
                            2.upto(n-1){|i|
                              thetaW += sfsArray[i] * i * (1.0 / i)
                            }
                            return(thetaW * (1.0 / weightSum))
                          end

                          ##Uses Achaz's "system" to calculate theta_pi based on sfs
                          # doesn't include singletons
                          def sfs2ThetaPiNoSingletons(sfsArray)
                            weightSum = 0.0
                            n = sfsArray.length
                            2.upto(n-1){|i|
                              weightSum += (n - i)
                            }
                            thetaW = 0.0
                            2.upto(n-1){|i|
                              thetaW += sfsArray[i] * i * (n- i)
                            }
                            return(thetaW * (1.0 / weightSum))
                          end
                          ##Uses Achaz's "system" to calculate theta_h based on sfs
                          # doesn't include singletons
                          def sfs2ThetaHNoSingletons(sfsArray)
                            weightSum = 0.0
                            n = sfsArray.length
                            2.upto(n-1){|i|
                              weightSum +=  i
                            }
                            thetaW = 0.0
                            2.upto(n-1){|i|
                              thetaW += sfsArray[i] * i *  i
                            }
                            return(thetaW * (1.0 / weightSum))
                          end

                          ##Uses Achaz's "system" to calculate Fay and Wu's H based on sfs
                          # doesn't include singletons
                          def sfs2FayWuHNoSingletons(sfsArray)
                            return(self.sfs2ThetaPiNoSingletons(sfsArray) - self.sfs2ThetaHNoSingletons(sfsArray) )
                          end

                          ##Uses Achaz's "system" to calculate theta_h based on sfs
                          # doesn't include singletons
                          def sfs2TajDNoSingletons(sfsArray)
                            return(self.sfs2ThetaPiNoSingletons(sfsArray) - self.sfs2ThetaWNoSingletons(sfsArray) )
                          end
                        end				

                        #Alignment.rb- a class for divergence statistics
                        #
                        #
                        class Alignment
                          attr_reader :sequenceMatrix1, :sequenceMatrix2, :outgroup, :codingRegions, :readingFrame, :geneticCode
                          attr_writer :sequenceMatrix1, :sequenceMatrix2, :outgroup, :codingRegions, :readingFrame, :geneticCode

                          def initialize(seqMat1, seqMat2, outG)
                            self.sequenceMatrix1 = seqMat1
                            self.sequenceMatrix2 = seqMat2
                            self.outgroup = outG
                            self.readingFrame = seqMat1.readingFrame		#inherits annotations from seqMat1
                            self.codingRegions = seqMat1.codingRegions
                            self.geneticCode = seqMat1.geneticCode
                          end


                          # manipulations
                          def length
                            return self.sequenceMatrix1.matrix[0].length
                          end

                          def alignSet(aSite)
                            a = Array.new
                            self.sequenceMatrix1.siteSet(aSite).each{ | x | a << x }
                            self.sequenceMatrix2.siteSet(aSite).each{ | x | a << x }
                            return a.uniq
                          end

                          def alignCodonSet(aCodonIndex)
                            a = Array.new
                            self.sequenceMatrix1.codonSet(aCodonIndex).each{ | x | a << x }
                            self.sequenceMatrix2.codonSet(aCodonIndex).each{ | x | a << x }
                            return a.uniq
                          end


                          def alignCodonSetPlus(aCodonIndex)
                            a = alignCodonSet(aCodonIndex)
                            a << self.outgroup.codonSet(aCodonIndex).first
                            return a.uniq
                          end

                          def alignAminoSet(aCodonIndex)
                            c = self.alignCodonSet(aCodonIndex)
                            c.reject!{ | x | x.include?("N")}
                            temp = Array.new
                            c.each{ | eachCodon |
                              temp << self.geneticCode[eachCodon]
                            }
                            return temp
                          end

                          def alignAminoSetPlus(aCodonIndex)
                            a = self.alignAminoSet(aCodonIndex)
                            a << self.outgroup.aminoSet(aCodonIndex).first
                            return a.uniq
                          end

                          def alignAminoSetPlusClean(aCodonIndex)
                            a = self.alignAminoSetPlus(aCodonIndex)
                            a.reject!{ | x | x == ("-") }
                            return a.uniq
                          end




                          #
                          # divergence
                          #
                          def countTransitions					#gets counts of differences according to basepair changes, only uses first sequence in each matrix
                            aHash = Hash.new
                            params = ["AT","AC","AG","CG","CT","GT","AA","CC","TT","GG"]
                            params.each{ | each | aHash[each] = 0 }
                            l = self.sequenceMatrix1.length - 1
                            0.upto(l){ | i | 
                              temp = Array.new
                              temp << self.sequenceMatrix1.matrix[0][i,1]
                              temp << self.sequenceMatrix2.matrix[0][i,1]
                              string = temp.sort.to_s
                              if params.include?(string)
                                aHash[string] += 1
                              end
                            }
                            params.each{ | each |
                              print each,"\t",aHash[each],"\n"
                            }
                          end


                          def fixedDiffSites(sampleSizeFilter)			#returns a dict with fixes and gaps/ambiguous
                            fixes = Array.new
                            gaps = Array.new
                            width = self.sequenceMatrix1.length
                            height1 = self.sequenceMatrix1.sampleSize
                            height2 = self.sequenceMatrix2.sampleSize
                            0.upto(width - 1){ | c |
                              sampSize = self.sequenceMatrix1.siteArrayClean(c).size
                              if sampSize > sampleSizeFilter
                                #use cleaned up data
                                testSet1 = self.sequenceMatrix1.siteSet(c).reject{| x | x == "N"}
                                testSet2 = self.sequenceMatrix2.siteSet(c).reject{| x | x == "N"}
                                if testSet1.empty? or testSet2.empty?
                                  gaps << c
                                else
                                  if testSet1.size == 1 and testSet2.size == 1 						#both sites fixed?
                                    if  (testSet1.include?("-") or testSet2.include?("-"))			#include gaps? if not...
                                      gaps << c	
                                    else
                                      if testSet1 != testSet2	  #are the sets equal? if not...
                                        fixes << c
                                      end
                                    end
                                  end
                                end
                              end
                            }
                            dict = Hash.new
                            dict["fixes"] = fixes
                            dict["gaps"] = gaps
                            return dict
                          end

                          def fixedDiffSitesWindow(start, fin)			#returns a dict with fixes and gaps/ambiguous
                            fixes = Array.new
                            gaps = Array.new
                            height1 = self.sequenceMatrix1.sampleSize
                            height2 = self.sequenceMatrix2.sampleSize
                            start.upto(fin - 1){ | c |
                              testSet1 = self.sequenceMatrix1.siteSet(c)
                              testSet2 = self.sequenceMatrix2.siteSet(c)
                              #clean up ambiguous sites
                              testSet1.delete("N")
                              testSet2.delete("N")
                              if testSet1.empty? or testSet2.empty?
                                gaps << c
                              else
                                if testSet1.size == 1 and testSet2.size == 1 						#both sites fixed?
                                  if  (testSet1.include?("-") or testSet2.include?("-"))			#include gaps? if not...
                                    gaps << c	
                                  else
                                    if testSet1 != testSet2										#are the sets equal? if not...
                                      fixes << c
                                    end
                                  end
                                end
                              end
                            }
                            dict = Hash.new
                            dict["fixes"] = fixes
                            dict["gaps"] = gaps
                            return dict
                          end

                          def fixedDiffsWindow(start, fin)										#returns number of fixed differences
                            fixes = 0
                            nCount = 0
                            start.upto(fin-1){ | c |
                              testSet1 = self.sequenceMatrix1.siteSet(c)
                              testSet2 = self.sequenceMatrix2.siteSet(c)
                              #clean up ambiguous sites
                              testSet1.delete("N")
                              testSet2.delete("N")
                              testSet1.delete("-")
                              testSet2.delete("-")
                              if testSet1.empty? or testSet2.empty?
                                nCount += 1
                              else
                                if testSet1.size == 1 and testSet2.size == 1 						#both sites fixed?
                                  # print testSet1[0],"\t",testSet2[0],"\n"
                                  if testSet1[0] != testSet2[0]										#are the sets equal? if not...
                                    fixes += 1
                                  end
                                end
                              end
                            }
                            if nCount > ((fin - start) * 0.9)
                              return "NA"
                            else
                              return fixes
                            end
                          end

                          def fixedDiffsWindowPer(start, fin)										#returns number of fixed differences
                            fixes = 0
                            nCount = 0
                            bi = 0
                            start.upto(fin - 1){ | c |
                              testSet1 = self.sequenceMatrix1.siteSet(c)
                              testSet2 = self.sequenceMatrix2.siteSet(c)
                              #clean up ambiguous sites
                              testSet1.delete("N")
                              testSet2.delete("N")
                              if testSet1.empty? or testSet2.empty?
                                nCount += 1
                              else
                                if testSet1.size == 1 and testSet2.size == 1 						#both sites fixed?
                                  if  (testSet1.include?("-") or testSet2.include?("-"))			#include gaps? if not...
                                  else
                                    if testSet1 != testSet2										#are the sets equal? if not...
                                      fixes += 1
                                    end
                                    bi += 1
                                  end
                                end
                              end
                            }
                            if nCount > ((fin - start) * 0.5)
                              return "NA"
                            else
                              return fixes.to_f / bi
                            end
                          end

                          def fixedDiffCodingDict
                            pos = Hash.new
                            fixDict = self.fixedDiffSites(1)   #currently only finds fixations with a sample size of 2 or greater in seqmat 1
                            fixSites = fixDict["fixes"]
                            codingSites = self.sequenceMatrix1.codingSites
                            fixSites = fixSites.select{ | each | codingSites.include?(each) }
                            codonPos = Array.new
                            fixSites.each{ | eachSite | codonPos << (codingSites.index(eachSite) / 3) }
                            pos["codonFixes"] = codonPos
                            pos["locusFixes"] = fixSites
                            pos["gaps"] = fixDict["gaps"]
                            return pos
                          end

                          def silReplFixations					#returns a hash of silent and replacement fixations
                            fixes = Hash.new
                            ["replacements","silents","missingData"].each{ | x | fixes[x] = Array.new }
                            dict = self.fixedDiffCodingDict
                            codingSites = self.sequenceMatrix1.codingSites
                            lastCodon = nil
                            dict["codonFixes"].each_with_index{ | site, index |
                              if site != lastCodon			#fixes the problem of multiple substitutions per codon
                                lastCodon = site
                                testSet = self.alignCodonSet(site)
                                if testSet.any?{ | x | x.include?("N") or x.include?("-") }		#clean up ambiguous data
                                  fixes["missingData"] << dict["locusFixes"][index]
                                else
                                  if self.sequenceMatrix1.codonSegSitesPositions(testSet).size == 1		#one fixed diff in the codon?
                                    amino1 = self.sequenceMatrix1.aminoSet(site)
                                    amino2 = self.sequenceMatrix2.aminoSet(site)
                                    if amino1.size < amino2.size
                                      if amino1.any?{ | aa | amino2.include?(aa) }
                                        fixes["silents"] << dict["locusFixes"].uniq[index]	
                                      else
                                        fixes["replacements"] << dict["locusFixes"].uniq[index]
                                      end
                                    else
                                      if amino2.any?{ | aa | amino1.include?(aa) }
                                        fixes["silents"] << dict["locusFixes"].uniq[index]
                                      else
                                        fixes["replacements"] << dict["locusFixes"].uniq[index]
                                      end
                                    end
                                  else				#more than one fixed difference in codon, use paths
                                    paths = Array.new
                                    set1 = self.sequenceMatrix1.codonSet(site)
                                    set2 = self.sequenceMatrix2.codonSet(site)
                                    set1.each{ | codon1 |
                                      set2.each{ | codon2 |
                                        print codon1,"\t",codon2,"\n"
                                        paths << self.sequenceMatrix1.lookupPath2Codons([codon1,codon2])
                                      }
                                    }
                                    #set up the path length variable
                                    short = 100
                                    shortPath = nil
                                    paths = paths.reject{ | x | x.nil? or x.include?("*") }

                                    #go through paths, find shortest with respect to Replacements
                                    paths.each{ | aPath |
                                      rLength = 0
                                      aPath.each{ | eachArray | if eachArray.include?("R")
                                        rLength += 1 
                                      end
                                    }
                                    if short > rLength
                                      short = rLength
                                      shortPath = aPath
                                    else
                                      if shortPath.size > aPath.size		# more replacements but less changes overall
                                        shortPath = aPath
                                      end
                                    end
                                  }

                                  shortPath.each{ | eachArray |
                                    if eachArray.first == "R"
                                      fixes["replacements"] << codingSites[ (site * 3) + eachArray.last.to_i]
                                    else
                                      if eachArray.first == "S"
                                        fixes["silents"] << codingSites[ (site * 3) + eachArray.last.to_i]
                                      end
                                    end
                                  }
                                end
                              end
                            end
                          }
                          fixes["missingData"] = fixes["missingData"] | dict["gaps"]
                          return fixes
                        end

                        def silReplFixationsGreedy	#returns a hash of silent and replacement fixations - greedy cleans out N's and Gaps
                          fixes = Hash.new
                          ["replacements","silents","missingData"].each{ | x | fixes[x] = Array.new }
                          dict = self.fixedDiffCodingDict
                          codingSites = self.sequenceMatrix1.codingSites
                          lastCodon = nil
                          dict["codonFixes"].each_with_index{ | site, index |
                            if site != lastCodon		#fixes the problem of multiple substitutions per codon
                              lastCodon = site
                              testSet = self.alignCodonSet(site)
                              testSet.reject!{ | x | x.include?("N") or x.include?("-") }		#clean up ambiguous data
                              if ! testSet.empty?
                                if self.sequenceMatrix1.codonSegSitesPositions(testSet).size == 1		#one fixed diff in the codon?
                                  amino1 = self.sequenceMatrix1.aminoSetClean(site)
                                  amino2 = self.sequenceMatrix2.aminoSetClean(site)
                                  if amino1.size < amino2.size
                                    if amino1.any?{ | aa | amino2.include?(aa) }
                                      fixes["silents"] << dict["locusFixes"].uniq[index]	
                                    else
                                      fixes["replacements"] << dict["locusFixes"].uniq[index]
                                    end
                                  else
                                    if amino2.any?{ | aa | amino1.include?(aa) }
                                      fixes["silents"] << dict["locusFixes"].uniq[index]
                                    else
                                      fixes["replacements"] << dict["locusFixes"].uniq[index]
                                    end
                                  end
                                else			#more than one fixed difference in codon, use paths
                                  paths = Array.new
                                  set1 = self.sequenceMatrix1.codonSetClean(site)
                                  set2 = self.sequenceMatrix2.codonSetClean(site)
                                  set1.each{ | codon1 |
                                    set2.each{ | codon2 |
                                      paths << self.sequenceMatrix1.lookupPath2Codons([codon1,codon2])
                                    }
                                  }
                                  #set up the path length variable
                                  short = 100
                                  shortPath = nil
                                  paths.reject!{ | x | x.nil? }
                                  paths.reject!{  | x | x.include?("*") }

                                  #go through paths, find shortest with respect to Replacements
                                  paths.each{ | aPath |
                                    rLength = 0
                                    aPath.each{ | eachArray | if eachArray.include?("R")
                                      rLength += 1 
                                    end
                                  }
                                  if short > rLength
                                    short = rLength
                                    shortPath = aPath
                                  else
                                    if shortPath.size > aPath.size		# more replacements but less changes overall
                                      shortPath = aPath
                                    end
                                  end
                                }
                                if shortPath != nil							#bug here!
                                  shortPath.each{ | eachArray |
                                    if eachArray.first == "R"
                                      fixes["replacements"] << codingSites[ (site * 3) + eachArray.last.to_i]
                                    else
                                      if eachArray.first == "S"
                                        fixes["silents"] << codingSites[ (site * 3) + eachArray.last.to_i]
                                      end
                                    end

                                  }
                                end
                              end
                            end
                          end
                        }
                        fixes["missingData"] = fixes["missingData"] | dict["gaps"]
                        return fixes
                      end
                      def polarizeAllFixations			#returns a hash of polarized fixations, designed for one outgroup
                        polarized = Hash.new
                        temp1 = Array.new
                        temp2 = Array.new
                        notPol = Array.new
                        sites = self.fixedDiffSites(1)									#use fixed differences between ingroups as starting material
                        sites["fixes"].each{ | eachSite |
                          os = self.outgroup.siteSet(eachSite)
                          if os != ["-"] 
                            set1 = self.sequenceMatrix1.siteSet(eachSite)			
                            set2 = self.sequenceMatrix2.siteSet(eachSite)			#skip those sites which have a gap in the outgroup
                            if os.any?{ | each | set1.include?(each) }				#outgroup includes set1?
                              if os.any?{ | each | set2.include?(each) }			#outgroup includes set2?
                                notPol << eachSite
                              else
                                temp2 << eachSite								# fixation along lineage 2
                              end
                            else
                              if os.any?{ | each | set2.include?(each) }								
                                temp1 << eachSite
                              end
                            end
                          else
                            notPol << eachSite
                          end
                        }
                        polarized["species1"] = temp1
                        polarized["species2"] = temp2
                        polarized["notPol"] = notPol
                        return polarized
                      end

                      def polarizeAllFixationsGreedy			#clean 'em up! ## returns a hash of polarized fixations, designed for one outgroup
                        polarized = Hash.new
                        temp1 = Array.new
                        temp2 = Array.new
                        notPol = Array.new
                        sites = self.fixedDiffSites(1)									#use fixed differences between ingroups as starting material
                        sites["fixes"].each{ | eachSite |
                          os = self.outgroup.siteSetClean(eachSite)
                          if os != ["-"] 
                            set1 = self.sequenceMatrix1.siteSetClean(eachSite)			
                            set2 = self.sequenceMatrix2.siteSetClean(eachSite)			#skip those sites which have a gap in the outgroup
                            if os.any?{ | each | set1.include?(each) }				#outgroup includes set1?
                              if os.any?{ | each | set2.include?(each) }			#outgroup includes set2?
                                notPol << eachSite
                              else
                                temp2 << eachSite								# fixation along lineage 2
                              end
                            else
                              if os.any?{ | each | set2.include?(each) }								
                                temp1 << eachSite
                              end
                            end
                          else
                            notPol << eachSite
                          end
                        }
                        polarized["species1"] = temp1
                        polarized["species2"] = temp2
                        polarized["notPol"] = notPol
                        return polarized
                      end
                      def polarizeCodingFixations
                        pol = Hash.new
                        polFixes = self.polarizeAllFixations											#use all polarized fixations
                        srHash = self.silReplFixations
                        codingSites = self.sequenceMatrix1.codingSites
                        pol["notPol"] = Array.new
                        ["species1","species2"].each{ | name |											#set up pol hash with an internal hash
                          aHash = Hash.new
                          ["silents","replacements","preferreds","unpreferreds","noChange"].each{ | x | aHash[x] = Array.new }
                          pol[name] = aHash
                        }

                        polFixes["species1"].each{ | eachFix |
                          if srHash["replacements"].include?(eachFix)									#silreplHash has site as repl
                            if self.alignAminoSetPlus(codingSites.index(eachFix) / 3).size < 3	#aminoSet can be polarized?
                              pol["species1"]["replacements"] << eachFix									
                            else
                              pol["notPol"] << eachFix
                            end
                          else
                            if srHash["silents"].include?(eachFix)
                              pol["species1"]["silents"] << eachFix								
                              #move to pref/unpref analysis
                              puSet1 = self.sequenceMatrix1.puSet(codingSites.index(eachFix) / 3)	
                              puSet2 = self.sequenceMatrix2.puSet(codingSites.index(eachFix) / 3)
                              puSetO = self.outgroup.puSet(codingSites.index(eachFix) / 3)
                              if ! (puSet1.empty? or puSet2.empty? or puSetO.empty?)
                                if puSet2 == puSet1
                                  pol["species1"]["noChange"] << eachFix
                                else
                                  if ! puSetO.include?(puSet1.first)						#is it derived?
                                    if puSet1.first == 1
                                      pol["species1"]["preferreds"] << eachFix
                                    else
                                      pol["species1"]["unpreferreds"] << eachFix
                                    end
                                  end
                                end
                              end
                            end
                          end
                        }
                        polFixes["species2"].each{ | eachFix |
                          if srHash["replacements"].include?(eachFix)									#silreplHash has site as repl
                            if self.alignAminoSetPlus( codingSites.index(eachFix) / 3 ).size < 3	#aminoSet can be polarized?
                              pol["species2"]["replacements"] << eachFix		
                            else
                              pol["notPol"] << eachFix
                            end
                          else
                            if srHash["silents"].include?(eachFix)
                              pol["species2"]["silents"] << eachFix								
                              #move to pref/unpref analysis
                              puSet1 = self.sequenceMatrix1.puSet(codingSites.index(eachFix) / 3)	
                              puSet2 = self.sequenceMatrix2.puSet(codingSites.index(eachFix) / 3)
                              puSetO = self.outgroup.puSet(codingSites.index(eachFix) / 3)
                              if ! (puSet1.empty? or puSet2.empty? or puSetO.empty?)
                                if puSet2 == puSet1
                                  pol["species2"]["noChange"] << eachFix
                                else
                                  if ! puSetO.include?(puSet2.first)						#is it derived?
                                    if puSet2.first == 1
                                      pol["species2"]["preferreds"] << eachFix
                                    else
                                      pol["species2"]["unpreferreds"] << eachFix
                                    end
                                  end
                                end
                              end
                            end
                          end
                        }
                        pol["notPol"] = pol["notPol"] | polFixes["notPol"]
                        return pol
                      end

                      def polarizeCodingFixationsGreedy			#clean 'em!
                        pol = Hash.new
                        polFixes = self.polarizeAllFixationsGreedy				#use all polarized fixations
                        srHash = self.silReplFixationsGreedy
                        codingSites = self.sequenceMatrix1.codingSites
                        pol["notPol"] = Array.new
                        ["species1","species2"].each{ | name |				#set up pol hash with an internal hash
                          aHash = Hash.new
                          ["silents","replacements","preferreds","unpreferreds","noChange"].each{ | x | aHash[x] = Array.new }
                          pol[name] = aHash
                        }

                        polFixes["species1"].each{ | eachFix |
                          if srHash["replacements"].include?(eachFix)									#silreplHash has site as repl
                            if self.alignAminoSetPlusClean(codingSites.index(eachFix) / 3).size < 3	#aminoSet can be polarized?
                              pol["species1"]["replacements"] << eachFix									
                            else
                              pol["notPol"] << eachFix
                            end
                          else
                            if srHash["silents"].include?(eachFix)
                              pol["species1"]["silents"] << eachFix								
                              #move to pref/unpref analysis
                              puSet1 = self.sequenceMatrix1.puSet(codingSites.index(eachFix) / 3)	
                              puSet2 = self.sequenceMatrix2.puSet(codingSites.index(eachFix) / 3)
                              puSetO = self.outgroup.puSet(codingSites.index(eachFix) / 3)
                              if ! (puSet1.empty? or puSet2.empty? or puSetO.empty?)
                                if puSet2 == puSet1
                                  pol["species1"]["noChange"] << eachFix
                                else
                                  if ! puSetO.include?(puSet1.first)						#is it derived?
                                    if puSet1.first == 1
                                      pol["species1"]["preferreds"] << eachFix
                                    else
                                      pol["species1"]["unpreferreds"] << eachFix
                                    end
                                  end
                                end
                              end
                            end
                          end
                        }
                        polFixes["species2"].each{ | eachFix |
                          if srHash["replacements"].include?(eachFix)									#silreplHash has site as repl
                            if self.alignAminoSetPlusClean( codingSites.index(eachFix) / 3 ).size < 3	#aminoSet can be polarized?
                              pol["species2"]["replacements"] << eachFix		
                            else
                              pol["notPol"] << eachFix
                            end
                          else
                            if srHash["silents"].include?(eachFix)
                              pol["species2"]["silents"] << eachFix								
                              #move to pref/unpref analysis
                              puSet1 = self.sequenceMatrix1.puSet(codingSites.index(eachFix) / 3)	
                              puSet2 = self.sequenceMatrix2.puSet(codingSites.index(eachFix) / 3)
                              puSetO = self.outgroup.puSet(codingSites.index(eachFix) / 3)
                              if ! (puSet1.empty? or puSet2.empty? or puSetO.empty?)
                                if puSet2 == puSet1
                                  pol["species2"]["noChange"] << eachFix
                                else
                                  if ! puSetO.include?(puSet2.first)						#is it derived?
                                    if puSet2.first == 1
                                      pol["species2"]["preferreds"] << eachFix
                                    else
                                      pol["species2"]["unpreferreds"] << eachFix
                                    end
                                  end
                                end
                              end
                            end
                          end
                        }
                        pol["notPol"] = pol["notPol"] | polFixes["notPol"]
                        return pol
                      end

                      #
                      #polymorphism
                      #

                      def silReplPolymorphismWithin		#returns a hash with silent and replacement polymorphisms that are specific to one of the seq mats
                        polys = Hash.new
                        ["replacements","silents","missingData"].each{|each | polys[each] = Array.new }
                        codingSites = self.sequenceMatrix1.codingSites
                        checkArray = Array.new
                        all = Array.new
                        hash1 = self.sequenceMatrix1.silReplMutations
                        hash2 = self.sequenceMatrix2.silReplMutations
                        hash1["silents"].each{ | each | 
                          all << each
                          if hash2["silents"].include?(each)
                            checkArray << each
                          end
                        }
                        hash2["silents"].each{ | each | all << each }
                        checkArray.uniq.each{ | x | all.delete(x) }
                        checkArray.uniq.each{ | site | 
                          codon1 = self.sequenceMatrix1.codonSet( codingSites.index(site) / 3 )
                          codon2 = self.sequenceMatrix2.codonSet( codingSites.index(site)/ 3 )
                          if codon1.size == codon2.size
                            if codon1.all?{ | each | codon2.include?(each) }
                              hash1["silents"].occurrencesOf(site).times{ all << site }
                            else
                              hash1["silents"].occurrencesOf(site).times{ all << site }
                              hash2["silents"].occurrencesOf(site).times{ all << site }
                            end
                          else
                            siteSet1 = self.sequenceMatrix1.siteSet(site)
                            siteSet2 = self.sequenceMatrix2.siteSet(site)
                            bothSet = Array.new
                            siteSet1.each{ | eachSite | bothSet << eachSite }
                            siteSet2.each{ | eachSite | bothSet << eachSite }
                            (bothSet.uniq.size - 1).times{ all << site }
                          end
                        }
                        all.each{ | eachSite | polys["silents"] << eachSite }

                        #add replacements

                        checkArray = Array.new
                        all = Array.new
                        hash1["replacements"].each{ | each | 
                          all << each
                          if hash2["replacements"].include?(each)
                            checkArray << each
                          end
                        }
                        hash2["replacements"].each{ | each | all << each }
                        checkArray.uniq.each{ | x | all.delete(x) }
                        checkArray.uniq.each{ | site | 
                          codon1 = self.sequenceMatrix1.codonSet( codingSites.index(site) / 3 )
                          codon2 = self.sequenceMatrix2.codonSet( codingSites.index(site) / 3 )
                          if codon1.size == codon2.size
                            if codon1.all?{ | each | codon2.include?(each) }
                              hash1["replacements"].occurrencesOf(site).times{ all << site }
                            else
                              hash1["replacements"].occurrencesOf(site).times{ all << site }
                              hash2["replacements"].occurrencesOf(site).times{ all << site }
                            end
                          else
                            siteSet1 = self.sequenceMatrix1.siteSet(site)
                            siteSet2 = self.sequenceMatrix2.siteSet(site)
                            bothSet = siteSet1 | siteSet2 
                            (bothSet.size - 1).times{ all << site }
                          end
                        }
                        all.each{ | eachSite | polys["replacements"] << eachSite }
                        polys["complex"] = hash1["complex"] | hash2["complex"]
                        return polys
                      end


                      def polarizeAllPolymorphisms
                        polys = Hash.new
                        ["species1","species2","notPol"].each{ | x | polys[x] = Array.new }
                        both = self.sequenceMatrix1.segSitesLocationsAllSites | self.sequenceMatrix2.segSitesLocationsAllSites
                        both.each{ | eachSite |
                          os = self.outgroup.siteSetClean(eachSite)
                          set1 = self.sequenceMatrix1.siteSetClean(eachSite)
                          set2 = self.sequenceMatrix2.siteSetClean(eachSite)
                          if set1.size == 1 or set2.size == 1				#one ingroup fixed?
                            if set1.size == 1					#set1 fixed? if true
                              if  set2.all?{ | each | os.include?(each) }	#outgroup has all set2 states? 
                                polys["notPol"]  << eachSite
                              else	
                                if set1.any?{ | each | os.include?(each) }		#outgroup has any set1 states? if true...
                                  polys["species2"] << eachSite 			#species2 polymorphism
                                else
                                  polys["notPol"]  << eachSite
                                end
                              end
                            else							#set2 fixed
                              if  set1.all?{ | each | os.include?(each) }
                                polys["notPol"]  << eachSite
                              else
                                if set2.any?{ | each | os.include?(each) }
                                  polys["species1"] << eachSite 
                                else
                                  polys["notPol"]  << eachSite
                                end
                              end
                            end
                          else								#both sets larger than 1 state
                            set1.each{ | x | 					# delete overlapping states
                              if set2.include?(x)
                                set1.delete(x)
                                set2.delete(x)
                              end
                            }
                            if os.any?{ | each | set1.include?(each) }			#outgroup state in set1? if true
                              if ! os.any?{ | each | set2.include?(each) }		#outgroup state in set2? if false
                                polys["species2"] << eachSite 			#species2 polymorphism
                              else 
                                polys["notPol"]  << eachSite			#outgroup has state in set1 and set2 -> nonPol
                              end
                            else
                              if  os.any?{ | each | set2.include?(each) }
                                polys["species1"] << eachSite 
                              end
                            end
                          end
                        }
                        return polys
                      end

                      def polarizeCodingPolymorphisms
                        pol = Hash.new
                        pol["notPol"] = Array.new
                        polPolys = self.polarizeAllPolymorphisms					#use all polarized polys
                        srHash = self.silReplPolymorphismWithin
                        codingSites = self.sequenceMatrix1.codingSites
                        ["species1","species2"].each{ | name |								#set up pol hash with an internal hash
                          aHash = Hash.new
                          ["silents","replacements","preferreds","unpreferreds","noChange"].each{ | x | aHash[x] = Array.new }
                          pol[name] = aHash
                        }
                        polPolys["species1"].each{ | eachPoly |
                          if srHash["replacements"].include?(eachPoly)									#silreplHash has site as repl
                            oAA = self.outgroup.aminoSet(codingSites.index(eachPoly) / 3)
                            inAA = self.sequenceMatrix2.aminoSet(codingSites.index(eachPoly) / 3)
                            if oAA.any?{ | each | inAA.include?(each) } 
                              pol["species1"]["replacements"] << eachPoly
                            else 
                              pol["notPol"] << eachPoly
                            end
                          else
                            if srHash["silents"].include?(eachPoly)
                              pol["species1"]["silents"] << eachPoly								
                              #move to pref/unpref analysis
                              puSet1 = self.sequenceMatrix1.puSet(codingSites.index(eachPoly) / 3)	
                              puSet2 = self.sequenceMatrix2.puSet(codingSites.index(eachPoly) / 3)
                              puSetO = self.outgroup.puSet(codingSites.index(eachPoly) / 3)
                              if ! (puSet1.empty? or puSet2.empty? or puSetO.empty?)
                                if puSet1.size == puSet2.size
                                  pol["species1"]["noChange"] << eachPoly
                                else
                                  puSet1 = puSet1.reject{ | x | puSet2.include?(x) } 
                                  if ! puSet1.empty?
                                    if puSet1.first == 1
                                      pol["species1"]["preferreds"] << eachPoly
                                    else
                                      pol["species1"]["unpreferreds"] << eachPoly
                                    end
                                  end
                                end
                              end
                            end
                          end
                        }
                        polPolys["species2"].each{ | eachPoly |
                          if srHash["replacements"].include?(eachPoly)						#silreplHash has site as repl
                            oAA = self.outgroup.aminoSet(codingSites.index(eachPoly) / 3)
                            inAA = self.sequenceMatrix1.aminoSet(codingSites.index(eachPoly) / 3)
                            if oAA.any?{ | each | inAA.include?(each) } 
                              pol["species2"]["replacements"] << eachPoly
                            else 
                              pol["notPol"] << eachPoly
                            end								
                          else
                            if srHash["silents"].include?(eachPoly)
                              pol["species2"]["silents"] << eachPoly								
                              #move to pref/unpref analysis
                              puSet1 = self.sequenceMatrix1.puSet(codingSites.index(eachPoly) / 3)	
                              puSet2 = self.sequenceMatrix2.puSet(codingSites.index(eachPoly) / 3)
                              puSetO = self.outgroup.puSet(codingSites.index(eachPoly) / 3)
                              if ! (puSet1.empty? or puSet2.empty? or puSetO.empty?)
                                if puSet1.size == puSet2.size
                                  pol["species2"]["noChange"] << eachPoly
                                else
                                  puSet2 = puSet2.reject{ | x | puSet1.include?(x) } 
                                  if ! puSet1.empty?
                                    if puSet2.first == 1
                                      pol["species2"]["preferreds"] << eachPoly
                                    else
                                      pol["species2"]["unpreferreds"] << eachPoly
                                    end
                                  end
                                end
                              end
                            end
                          end
                        }
                        pol["notPol"] = pol["notPol"] | polPolys["notPol"]
                        return pol
                      end

                      def fst
                        n = self.sequenceMatrix1.sampleSize + self.sequenceMatrix2.sampleSize
                        hb = self.averagePairwiseDist 
                        print hb,"\t"
                        hw = self.sequenceMatrix1.averagePairwiseDifferencesArray  + self.sequenceMatrix2.averagePairwiseDifferencesArray  

                        print hw.sampleMean,"\n"
                        if hw == "NA" or hb == "NA"
                          return("NA")
                        else
                          return(1.0 - (hw.sampleMean/hb))
                        end
                      end


                      def averagePairwiseDist
                        ds = Array.new      
                        d = 0
                        #go through all pairwise combos
                        self.sequenceMatrix1.matrix.each_index{ | i |
                          self.sequenceMatrix2.matrix.each_index{ | j |
                            #print i,"\t",j,"\n"
                            tempSeq1 = self.sequenceMatrix1.matrix[i]
                            tempSeq2 = self.sequenceMatrix2.matrix[j]
                            ds << Alignment.stringComp(tempSeq1,tempSeq2)
                          }
                        }
                        return ds.sampleMean
                      end

                      def averagePairwiseJCDist
                        ds = Array.new

                        #go through all pairwise combos
                        self.sequenceMatrix1.matrix.each_index{ | i |
                          self.sequenceMatrix2.matrix.each_index{ | j |
                            tempSeq1 = self.sequenceMatrix1.returnSequences([i])
                            tempSeq2 = self.sequenceMatrix2.returnSequences([j])
                            tempAlign = Alignment.new(tempSeq1,tempSeq2,nil)
                            d = tempAlign.jcDiv
                            ds << d

                          }
                        }
                        return ds.sampleMean
                      end

                      def averagePairwiseNGDist     #returns an array (dn,ds), this is algorithm 1
                        dns = Array.new
                        dss = Array.new

                        #go through all the pairwise combos
                        self.sequenceMatrix1.matrix.each_index{ | i |
                          self.sequenceMatrix2.matrix.each_index{ | j |
                            dists = silReplDifsAllPaths(i,j)
                            if ! dists.include?(nil)
                              dns << (-3.0/4.0) * Math.log(1.0 - ((4.0/3.0)*dists[0]))
                              dss << (-3.0/4.0) * Math.log(1.0 - ((4.0/3.0)*dists[1]))
                            end
                          }
                        }
                        return [dns.sampleMean, dss.sampleMean]
                      end
                      def averageWeightedPairwiseNGDist     #same as above but weights final average by sil/repl sites
                        dns = Array.new
                        dss = Array.new

                        #make arrays to store average sites for comparisons
                        ssComps = Array.new
                        rsComps = Array.new
                        #go through all the pairwise combos
                        self.sequenceMatrix1.matrix.each_index{ | i |
                          self.sequenceMatrix2.matrix.each_index{ | j |
                            dists = silReplDifsAllPaths(i,j)
                            if ! dists.include?(nil)
                              sites = silReplSitesInComp(i,j)
                              rsComps << sites[0]
                              ssComps << sites[1]
                              dns << (-3.0/4.0) * Math.log(1.0 - ((4.0/3.0)*dists[0]))
                              dss << (-3.0/4.0) * Math.log(1.0 - ((4.0/3.0)*dists[1]))
                            end
                          }
                        }
                        #calc weights
                        ssSum = ssComps.sum
                        rsSum = rsComps.sum
                        ssComps.each_index{ | i |
                          dss[i] = dss[i] * (ssComps[i] / ssSum)
                          dns[i] = dns[i] * (rsComps[i] / rsSum)
                        }

                        return [dns.sum, dss.sum]
                      end
                      def averagePairwiseNGDist2     #returns an array (dn,ds), this is algorithm 2
                        dns = Array.new
                        dss = Array.new

                        #go through all the pairwise combos
                        self.sequenceMatrix1.matrix.each_index{ | i |
                          self.sequenceMatrix2.matrix.each_index{ | j |
                            dists = silReplDifsNG(i,j)
                            dns << (-3.0/4.0) * Math.log(1.0 - ((4.0/3.0)*dists[0]))
                            dss << (-3.0/4.0) * Math.log(1.0 - ((4.0/3.0)*dists[1]))
                          }
                        }
                        return [dns.sampleMean, dss.sampleMean]
                      end
                      def averageWeightedPairwiseNGDist2     #same as above but weights final average by sil/repl sites
                        dns = Array.new
                        dss = Array.new
                        #get silent and replacement sites for each species
                        ss1 = self.sequenceMatrix1.silentSites
                        rs1 = self.sequenceMatrix1.replacementSites
                        ss2 = self.sequenceMatrix2.silentSites
                        rs2 = self.sequenceMatrix2.replacementSites
                        #make arrays to store average sites for comparisons
                        ssComps = Array.new
                        rsComps = Array.new
                        #go through all the pairwise combos
                        self.sequenceMatrix1.matrix.each_index{ | i |
                          self.sequenceMatrix2.matrix.each_index{ | j |
                            numerator = silReplDifsNG(i,j)
                            denom = [rs1[i],rs2[j]].sampleMean
                            rsComps << denom
                            dns << (-3.0/4.0) * Math.log(1.0 - ((4.0/3.0)*(numerator[0].to_f / denom)))
                            denom = [ss1[i],ss2[j]].sampleMean
                            ssComps << denom
                            dss << (-3.0/4.0) * Math.log(1.0 - ((4.0/3.0)*(numerator[1].to_f / denom)))

                          }
                        }
                        #calc weights
                        ssSum = ssComps.sum
                        rsSum = rsComps.sum
                        ssComps.each_index{ | i |
                          dss[i] = dss[i] * (ssComps[i] / ssSum)
                          dns[i] = dns[i] * (rsComps[i] / rsSum)
                        }

                        return [dns.sum, dss.sum]
                      end


                      def silReplDifsNG(specAIndex, specBIndex)  #this returns an array (nonSyn, Syn) as estimated by Nei & Gojobori alg. 2

                        codons1 = self.sequenceMatrix1.codons[specAIndex]
                        codons2 = self.sequenceMatrix2.codons[specBIndex]

                        #initialize counts
                        ms = [0,0,0]
                        mm = [0,0,0]
                        s = [0,0,0]

                        #go through codons, count up number of diffs in each codon and their positions
                        codons1.each_index{ | i |
                          set = Array.new 
                          set << codons1[i]
                          set << codons2[i]
                          #get rid of any ambiguous codons
                          set.reject!{ | x | x.include?("N") or x.include?("X") }
                          ssPos = self.sequenceMatrix1.codonSegSitesPositions(set)
                          #one change?
                          if ssPos.size == 1
                            ms[ssPos.first] += 1
                            #is it silent?
                            if self.sequenceMatrix1.aminoCodonSet(set).size == 1
                              s[ssPos.first] += 1
                            end
                          else
                            #two or three changes in codon?
                            if ssPos.size > 1
                              ssPos.each{ | j| 
                                mm[j] += 1
                              }
                            end
                          end
                        }
                        #get prob of synonymous per position
                        if s[0] == 0
                          probSyn1 = 0.0
                        else
                          probSyn1 = s[0].to_f / ms[0]
                        end

                        if s[2] == 0
                          probSyn3 = 0.0
                        else
                          probSyn3 = s[2].to_f / ms[2]
                        end
                        probNonSyn1 = 1.0 - probSyn1
                        probNonSyn3 = 1.0 - probSyn3

                        #estimate numbers
                        sd = (probSyn1 * (ms[0] + mm[0])) + (probSyn3 * (ms[2] + mm[2]))
                        nd = (probNonSyn1 * (ms[0] + mm[0])) +  (ms[1] + mm[1]) + (probNonSyn3 * (ms[2] + mm[2])) #assumes all second positions are nonsyn
                        return [nd,sd]
                      end

                      def silReplDifsAllPaths(specAIndex, specBIndex)  #this returns an array (nonSyn, Syn) as estimated by Nei & Gojobori alg. 1

                        codons1 = self.sequenceMatrix1.codons[specAIndex]
                        codons2 = self.sequenceMatrix2.codons[specBIndex]

                        #initialize arrays
                        nsCount = 0
                        sCount = 0
                        nsSites = 0
                        sSites = 0
                        #go through codons, count up number of diffs in each codon and their positions
                        codons1.each_index{ | i |
                          set = Array.new 
                          set << codons1[i]
                          set << codons2[i]
                          #get rid of any ambiguous codons
                          set.reject!{ | x | x.include?("N") or x.include?("X") or self.geneticCode[x] == "*" }
                          if set.size > 1
                            #tally up sites along the way
                            nsSites += self.sequenceMatrix1.silentSiteDict[set[0]][1]
                            nsSites += self.sequenceMatrix1.silentSiteDict[set[1]][1]
                            sSites += self.sequenceMatrix1.silentSiteDict[set[0]][0]
                            sSites += self.sequenceMatrix1.silentSiteDict[set[1]][0]
                          end
                          #more than one codon state?
                          if set.uniq.size > 1
                            #get allPaths between codons
                            paths = self.sequenceMatrix1.allPaths2Codons(set)
                            tempRCount = 0
                            tempSCount = 0
                            #clean out paths which go through stop codons
                            paths.reject!{ | x | x[0] == "*"}
                            #go through paths, count stuff
                            paths.each{ | aPath |
                              #go through path steps, count em
                              aPath.each{ | aStep |
                                if aStep[0] == "R"
                                  tempRCount += 1
                                else
                                  tempSCount += 1
                                end
                              }
                            }
                            nsCount += tempRCount.to_f / paths.size
                            sCount += tempSCount.to_f/ paths.size
                          end
                        }

                        if (nsSites == 0 and sSites == 0)
                          return [nil,nil]
                        else
                          return [nsCount / (nsSites.to_f / 2), sCount / (sSites.to_f / 2)]    
                        end
                      end
                      def silReplSitesInComp(specAIndex, specBIndex)  #this returns an array (nonSyn, Syn) of sites used in comparison

                        codons1 = self.sequenceMatrix1.codons[specAIndex]
                        codons2 = self.sequenceMatrix2.codons[specBIndex]

                        #initialize arrays
                        nsCount = 0
                        sCount = 0
                        nsSites = 0
                        sSites = 0
                        #go through codons, count up number of diffs in each codon and their positions
                        codons1.each_index{ | i |
                          set = Array.new 
                          set << codons1[i]
                          set << codons2[i]
                          #get rid of any ambiguous codons
                          set.reject!{ | x | x.include?("N") or x.include?("X") or self.sequenceMatrix1.geneticCode[x] == "*"}

                          #more than one codon state?
                          if set.size > 1
                            #tally up sites along the way
                            nsSites += self.sequenceMatrix1.silentSiteDict[set[0]][1]
                            nsSites += self.sequenceMatrix1.silentSiteDict[set[1]][1]
                            sSites += self.sequenceMatrix1.silentSiteDict[set[0]][0]
                            sSites += self.sequenceMatrix1.silentSiteDict[set[1]][0]
                          end
                        }
                        if (nsSites == 0 and sSites == 0)
                          return [nil,nil]
                        else
                          return [nsSites.to_f / 2, sSites.to_f / 2]    
                        end
                      end                        


                      #statistics


                      def mkTest				#returns an array in the order aaFix, aaPoly, silFix, silPoly
                        fix = self.silReplFixations
                        poly = self.silReplPolymorphismWithin
                        return [fix["replacements"].size,poly["replacements"].size,fix["silents"].size,poly["silents"].size]
                      end

                      def mkTestGreedy				#returns an array in the order aaFix, aaPoly, silFix, silPoly
                        fix = self.silReplFixationsGreedy
                        poly = self.silReplPolymorphismWithin
                        return [fix["replacements"].size,poly["replacements"].size,fix["silents"].size,poly["silents"].size]
                      end

                      def polMKTest			#returns hash
                        fix = self.polarizeCodingFixations
                        poly = self.polarizeCodingPolymorphisms
                        dict = Hash.new
                        ["species1","species2"].each{ | x |
                          temp = Hash.new
                          temp["aaFix"] = fix[x]["replacements"].size
                          temp["silFix"] = fix[x]["silents"].size
                          temp["prefFix"] = fix[x]["preferreds"].size
                          temp["unFix"] = fix[x]["unpreferreds"].size
                          temp["ncFix"] = fix[x]["noChange"].size
                          temp["aaPoly"] = poly[x]["replacements"].size
                          temp["silPoly"] = poly[x]["silents"].size
                          temp["prefPoly"] = poly[x]["preferreds"].size
                          temp["unPoly"] = poly[x]["unpreferreds"].size
                          temp["ncPoly"] = poly[x]["noChange"].size
                          dict[x] = temp
                        }
                        return dict
                      end

                      def polMKTestGreedy			#returns hash
                        fix = self.polarizeCodingFixationsGreedy
                        poly = self.polarizeCodingPolymorphisms
                        dict = Hash.new
                        ["species1","species2"].each{ | x |
                          temp = Hash.new
                          temp["aaFix"] = fix[x]["replacements"].size
                          temp["silFix"] = fix[x]["silents"].size
                          temp["prefFix"] = fix[x]["preferreds"].size
                          temp["unFix"] = fix[x]["unpreferreds"].size
                          temp["ncFix"] = fix[x]["noChange"].size
                          temp["aaPoly"] = poly[x]["replacements"].size
                          temp["silPoly"] = poly[x]["silents"].size
                          temp["prefPoly"] = poly[x]["preferreds"].size
                          temp["unPoly"] = poly[x]["unpreferreds"].size
                          temp["ncPoly"] = poly[x]["noChange"].size
                          dict[x] = temp
                        }
                        return dict
                      end
                      def polMKTestArray 		#returns 2 dimensional mkArray: aaFix, aaPoly, silFix, silPoly
                        polMK = polMKTest
                        out = Array.new
                        out << [polMK["species1"]["aaFix"],polMK["species1"]["aaPoly"],polMK["species1"]["silFix"],polMK["species1"]["silPoly"]]
                        out << [polMK["species2"]["aaFix"],polMK["species2"]["aaPoly"],polMK["species2"]["silFix"],polMK["species2"]["silPoly"]]
                        return out
                      end
                      def polMKTestArrayGreedy 		#returns 2 dimensional mkArray: aaFix, aaPoly, silFix, silPoly
                        polMK = polMKTestGreedy
                        out = Array.new
                        out << [polMK["species1"]["aaFix"],polMK["species1"]["aaPoly"],polMK["species1"]["silFix"],polMK["species1"]["silPoly"]]
                        out << [polMK["species2"]["aaFix"],polMK["species2"]["aaPoly"],polMK["species2"]["silFix"],polMK["species2"]["silPoly"]]
                        return out
                      end
                      def polPUTestArray 		#returns 2 dimensional mkArray: aaFix, aaPoly, silFix, silPoly
                        polMK = polMKTest
                        out = Array.new
                        out << [polMK["species1"]["prefFix"],polMK["species1"]["prefPoly"],polMK["species1"]["unFix"],polMK["species1"]["unPoly"],polMK["species1"]["ncFix"],polMK["species1"]["ncPoly"]]
                        out << [polMK["species2"]["prefFix"],polMK["species2"]["prefPoly"],polMK["species2"]["unFix"],polMK["species2"]["unPoly"],polMK["species2"]["ncFix"],polMK["species2"]["ncPoly"]]
                        return out
                      end

                      #
                      #	
                      #misc
                      #
                      #
                      def padAlignment		#evens sequence length with N's from the back
                        max = 0
                        self.sequenceMatrix1.matrix.each{ | eachString |
                          if eachString.length > max
                            max = eachString.length
                          end
                        }
                        self.sequenceMatrix2.matrix.each{ | eachString |
                          if eachString.length > max
                            max = eachString.length
                          end
                        }
                        self.sequenceMatrix1.matrix.each{ | eachString |
                          if eachString.length < max
                            l = max - eachString.length
                            l.times{ | i | eachString << "N" }
                          end
                        }
                        self.sequenceMatrix2.matrix.each{ | eachString |
                          if eachString.length < max
                            l = max - eachString.length
                            l.times{ | i | eachString << "N" }
                          end
                        }	
                      end
                      def slidingWindowFixedDiffs(windowSize, offset)    # [site, fd]
                        i = 0
                        l = self.length - windowSize
                        oc = Array.new
                        while i <= l
                          oc << [i, self.fixedDiffsWindow(i, i + windowSize)]
                          i += offset
                        end
                        oc.reject!{ | each | each[1] == "NA"}
                        oc.each{ | each | print each[0],"\t",each[1].to_f / windowSize,"\n"}
                      end

                      def slidingWindowJCDiv(windowSize, offset)    # only use with two sequences!!! [site, fd]
                        i = 0
                        l = self.length - windowSize
                        oc = Array.new
                        while i <= l
                          d =  self.fixedDiffsWindowPer(i, i + windowSize)
                          if d != "NA"
                            oc << [i, (-3.0/4) * Math.log(1 - ((4.0 *d )/3))]
                          end
                          i += offset
                        end
                        oc.each{ | each | print each[0],"\t",each[1],"\n"}
                      end

                      def jcDiv
                        d =  self.fixedDiffsWindowPer(0, self.length - 1)
                        if d != "NA"
                          return  (-3.0/4) * Math.log(1 - ((4.0 *d )/3))
                        else
                          return "NA"
                        end
                      end

                      #class methods
                      def Alignment.stringComp(string1,string2)
                        count = 0
                        0.upto(string1.length-1){|i|
                          if string1[i,1] != "N" and string2[i,1] != "N" and string1[i,1] != "-" and string2[i,1] != "-"
                            if string1[i,1] != string2[i,1]
                              count += 1
                            end
                          end
                        }
                        return(count)
                      end

                      def Alignment.gTest(anArray)		#takes an array of length 4, [a,c,b,d]
                        cellSum = (anArray[0] * Math.log(anArray[0])) + (anArray[1] * Math.log(anArray[1])) + (anArray[2] * Math.log(anArray[2])) + (anArray[3] * Math.log(anArray[3]))   
                        r1 = anArray[0] + anArray[2]
                        r2 = anArray[1] + anArray[3]
                        c1 = anArray[0] + anArray[1]
                        c2 = anArray[2] + anArray[3]
                        n = anArray.sum.to_f
                        nln = n * Math.log(n)
                        rcSum = (r1 * Math.log(r1)) + (r2 * Math.log(r2)) + (c1 * Math.log(c1)) + (c2 * Math.log(c2))
                        g = 2 * (cellSum - rcSum + nln)
                        q = 1 + ((((n/r1) + (n/r2) - 1) * ((n/c1)+(n/c2) - 1)) / (6 * n))
                        return g/q
                      end

                      def Alignment.fishersExactTest(anArray)		#takes an array of length 4, [a,c,b,d]
                        tmpArray= Array.new
                        anArray.each{ | x | tmpArray << x}
                        n = anArray.sum.to_i
                        q1 = (anArray[0] + anArray[2]).logFactorial + (anArray[1] + anArray[3]).logFactorial + (anArray[0] + anArray[1]).logFactorial + (anArray[2] + anArray[3]).logFactorial - n.logFactorial
                        q2 = anArray[0].logFactorial + anArray[1].logFactorial + anArray[2].logFactorial + anArray[3].logFactorial
                        pTail1 = 10 ** (q1 - q2)
                        origP1 = 10 ** (q1 - q2)
                        while (! anArray.include?(0))
                          if (anArray[0] * anArray[3]) - (anArray[1] * anArray[2]) < 0
                            anArray[0] -= 1
                            anArray[3] -= 1
                            anArray[1] += 1
                            anArray[2] += 1
                          else
                            anArray[0] += 1
                            anArray[3] += 1
                            anArray[1] -= 1
                            anArray[2] -= 1
                          end
                          q2 = anArray[0].logFactorial + anArray[1].logFactorial + anArray[2].logFactorial + anArray[3].logFactorial
                          pTail1 += 10 ** (q1 - q2)

                        end
                        #now add up second tail
                        if (tmpArray[0] * tmpArray[3]) - (tmpArray[1] * tmpArray[2]) < 0
                          adj = [anArray[1],anArray[2]].min
                          tmpArray[1] -= adj
                          tmpArray[2] -= adj
                          tmpArray[0] += adj
                          tmpArray[3] += adj
                        else
                          adj = [anArray[0],anArray[3]].min
                          tmpArray[1] += adj
                          tmpArray[2] += adj
                          tmpArray[0] -= adj
                          tmpArray[3] -= adj
                        end
                        q1 = (tmpArray[0] + tmpArray[2]).logFactorial + (tmpArray[1] + tmpArray[3]).logFactorial + (tmpArray[0] + tmpArray[1]).logFactorial + (tmpArray[2] + tmpArray[3]).logFactorial - n.logFactorial
                        q2 = tmpArray[0].logFactorial + tmpArray[1].logFactorial + tmpArray[2].logFactorial + tmpArray[3].logFactorial
                        pTail2 = 10 ** (q1 - q2)
                        origP2 = 10 ** (q1 - q2)
                        while (origP2 < origP1)
                          if (tmpArray[0] * tmpArray[3]) - (tmpArray[1] * tmpArray[2]) < 0
                            tmpArray[0] += 1
                            tmpArray[3] += 1
                            tmpArray[1] -= 1
                            tmpArray[2] -= 1
                          else
                            tmpArray[0] -= 1
                            tmpArray[3] -= 1
                            tmpArray[1] += 1
                            tmpArray[2] += 1
                          end
                          q2 = tmpArray[0].logFactorial + tmpArray[1].logFactorial + tmpArray[2].logFactorial + tmpArray[3].logFactorial
                          pTail2 += 10 ** (q1 - q2)
                          origP2 = 10 ** (q1 - q2)
                        end
                        pTail2 -= origP2
                        return pTail1 + pTail2
                      end
                    end	

###############################
########
#########           Main Loop for MK test below
#####
########
######
############################################

if ARGV.empty?
  print "mkTest.rb ingroup.fa outgroup.fa\n"
  print "\toptions:\n\t\t-p outgroup2.fa (polarized MK test)\n"
  exit
end

ing = SequenceMatrix.new.initializeFromFasta(ARGV[0])
outg = SequenceMatrix.new.initializeFromFasta(ARGV[1])

if ARGV.include?("-p")
  #Polarized MK test
  outg2 = SequenceMatrix.new.initializeFromFasta(ARGV[ARGV.index("-p")+1])
  ing.asCodingSequence(nil,nil)
  outg.asCodingSequence(nil,nil)
  ing.asCodingSequence(nil,nil)
  outg2.asCodingSequence(nil,nil)
  align = Alignment.new(ing,outg,outg2)
  array = align.polMKTestArrayGreedy
  print "popn1 aaFix\tpopn1 aaPoly\tpopn1 silFix\tpopn1 silPoly\tpopn2 aaFix\tpopn2 aaPoly\tpopn2 silFix\tpopn2 silPoly\n"
  print array[0][0],"\t",array[0][1],"\t",array[0][2],"\t",array[0][3],"\t"
	print array[1][0],"\t",array[1][1],"\t",array[1][2],"\t",array[1][3],"\n"


else
  #normal MK test
  print "aaFix\taaPoly\tsilFix\tsilPoly\n"
  ing.asCodingSequence(nil,nil)
  outg.asCodingSequence(nil,nil)
  align = Alignment.new(ing,outg,nil)

  #print align.sequenceMatrix1.matrix[0],"\n"
  array = align.mkTestGreedy
  print array[0],"\t",array[1],"\t",array[2],"\t",array[3],"\n"

end


